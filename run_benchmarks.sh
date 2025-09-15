#!/bin/bash

# AmpliconHunter Comprehensive Benchmarking Suite - All Versions
# Usage: ./run_benchmarks_all_versions.sh <results_directory> [test_number|all]
#   test_number: 1, 2, 3, 4, 5, or comma-separated list (e.g., "1,3,5")
#   default: all

set -e

# Parse arguments
if [ $# -lt 1 ] || [ $# -gt 2 ]; then
    echo "Usage: $0 <results_directory> [test_number|all]"
    echo "  test_number: 1, 2, 3, 4, 5, 6, or comma-separated list (e.g., '1,3,5')"
    echo "  default: all"
    echo ""
    echo "Tests:"
    echo "  1: Input Size Scaling"
    echo "  2: Primer N-base Experiment"
    echo "  3: Mismatch Variation"
    echo "  4: Thread Scaling"
    echo "  5: Hot vs Cold Cache"
    echo "  6: Primer Pair Comparison"
    echo ""
    echo "Examples:"
    echo "  $0 results         # Run all tests"
    echo "  $0 results all     # Run all tests"
    echo "  $0 results 5       # Run only test 5"
    echo "  $0 results 1,4,5   # Run tests 1, 4, and 5"
    exit 1
fi

RESULTS_DIR="$1"
TESTS_TO_RUN="${2:-all}"

# Parse which tests to run
RUN_TEST1=false
RUN_TEST2=false
RUN_TEST3=false
RUN_TEST4=false
RUN_TEST5=false
RUN_TEST6=false

if [ "$TESTS_TO_RUN" = "all" ]; then
    RUN_TEST1=true
    RUN_TEST2=true
    RUN_TEST3=true
    RUN_TEST4=true
    RUN_TEST5=true
    RUN_TEST6=true
else
    # Parse comma-separated list
    IFS=',' read -ra TESTS <<< "$TESTS_TO_RUN"
    for test in "${TESTS[@]}"; do
        case "$test" in
            1) RUN_TEST1=true ;;
            2) RUN_TEST2=true ;;
            3) RUN_TEST3=true ;;
            4) RUN_TEST4=true ;;
            5) RUN_TEST5=true ;;
            6) RUN_TEST6=true ;;
            *) echo "Invalid test number: $test"; exit 1 ;;
        esac
    done
fi

echo "=========================================="
echo "AmpliconHunter Comprehensive Benchmarking Suite"
echo "=========================================="
echo "Results directory: $RESULTS_DIR"
echo "Tests to run:"
[ "$RUN_TEST1" = true ] && echo "  - Test 1: Input Size Scaling"
[ "$RUN_TEST2" = true ] && echo "  - Test 2: Primer N-base Experiment"
[ "$RUN_TEST3" = true ] && echo "  - Test 3: Mismatch Variation"
[ "$RUN_TEST4" = true ] && echo "  - Test 4: Thread Scaling"
[ "$RUN_TEST5" = true ] && echo "  - Test 5: Hot vs Cold Cache"
[ "$RUN_TEST6" = true ] && echo "  - Test 6: Primer Pair Comparison"
echo "=========================================="
echo ""

BASE_DIR="$(pwd)"
DATA_DIR="${BASE_DIR}/data"
COMPRESSED_DIR="${DATA_DIR}/compressed"
FASTA_DIR="${DATA_DIR}/fasta"
PRIMERS="${DATA_DIR}/primers.txt"

# All versions to test
# Comment the line below to not test all versions:
VERSIONS=("AHv1" "AHv2" "AHv3" "AHv4" "AHv5")
# And uncomment this to test only specific versions:
# VERSIONS=("AHv5")

# Create results directory structure for all versions
for test_dir in test1 test2 test3 test4 test5 test6; do
    for version in "${VERSIONS[@]}"; do
        mkdir -p "${RESULTS_DIR}/${test_dir}/${version}"
    done
done

# Settings for all tests
THREADS=$(nproc)  # Use actual CPU cores
echo "Using ${THREADS} threads (based on CPU cores)"

CLAMP=3
MIN_LENGTH=1000
MAX_LENGTH=2000
FB_LEN=16
RB_LEN=16
REPEATS=3
TIMEOUT_SECONDS=3600  # 1 hour timeout
COOLDOWN_SECONDS=30   # Cooldown between runs to avoid thermal issues

# Set CPU governor to performance mode (requires sudo)
set_performance_governor() {
    if sudo -n true 2>/dev/null; then
        echo "Setting CPU governor to performance mode..."
        sudo bash -c 'for c in /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor; do echo performance > "$c" 2>/dev/null || true; done'
        echo "Current CPU frequencies:"
        awk '/cpu MHz/ {print}' /proc/cpuinfo | head -5
    else
        echo "WARNING: Cannot set CPU governor (no passwordless sudo). Results may vary due to frequency scaling."
    fi
}

# Proper cache clearing with sync
clear_cache() {
    echo "Clearing cache..."
    # First sync to flush dirty pages
    sync
    
    # Clear any leftover temp files
    rm -rf /tmp/amplicon_hunter_* 2>/dev/null || true
    
    # Check if we can use sudo without password
    if sudo -n true 2>/dev/null; then
        echo 3 | sudo tee /proc/sys/vm/drop_caches > /dev/null
        echo "Cache cleared successfully"
    else
        echo "WARNING: Cannot clear cache (no passwordless sudo). Results may be affected by caching."
    fi
    
    # Show cache status
    echo "Memory status:"
    grep -E 'Cached|Buffers|Dirty' /proc/meminfo
    
    sleep 2
}

# Function to convert time format M:SS.ss to seconds
convert_time_to_seconds() {
    local time_str="$1"
    if [[ "$time_str" == *":"* ]]; then
        # Format is M:SS.ss or H:MM:SS
        local parts=(${time_str//:/ })
        if [ ${#parts[@]} -eq 3 ]; then
            # H:MM:SS
            echo "scale=2; ${parts[0]} * 3600 + ${parts[1]} * 60 + ${parts[2]}" | bc
        else
            # M:SS.ss
            echo "scale=2; ${parts[0]} * 60 + ${parts[1]}" | bc
        fi
    else
        # Already in seconds
        echo "$time_str"
    fi
}

# Function to compress genomes for versions that need it
compress_genomes() {
    local version="$1"
    local genome_dir="$2"
    local output_dir="$3"
    local thread_count="${4:-$THREADS}"  # Optional thread count parameter
    
    # Only AHv2, AHv3, AHv4, and AHv5 need compression
    if [[ "${version}" != "AHv2" && "${version}" != "AHv3" && "${version}" != "AHv4" && "${version}" != "AHv5" ]]; then
        echo "Skipping compression for ${version}"
        return 0
    fi
    
    # AHv3-5 share the same compression
    if [ "${version}" == "AHv3" ]; then
        # Use AHv4's compressed data
        output_dir="${output_dir/AHv3/AHv4}"
    fi
    if [ "${version}" == "AHv5" ]; then
        # Use AHv4's compressed data
        output_dir="${output_dir/AHv5/AHv4}"
    fi
    
    if [ -d "${output_dir}" ] && [ "$(ls -A ${output_dir} 2>/dev/null)" ]; then
        echo "Compressed data already exists in ${output_dir}, skipping compression"
        return 0
    fi
    
    mkdir -p "${output_dir}"
    
    # Create compression timing directory
    local COMPRESS_LOG_DIR="${RESULTS_DIR}/compression_timing"
    mkdir -p "${COMPRESS_LOG_DIR}"
    
    # Generate log file name based on dataset size
    local DATASET_NAME=$(basename "${genome_dir}")
    if [ "${thread_count}" != "${THREADS}" ]; then
        DATASET_NAME="${DATASET_NAME}_t${thread_count}"
    fi
    local LOG_FILE="${COMPRESS_LOG_DIR}/${version}_${DATASET_NAME}_compression.log"
    
    # Clear cache before compression
    clear_cache
    
    # Determine which binary to use for compression
    local COMPRESS_BINARY
    case "${version}" in
        "AHv2")
            COMPRESS_BINARY="./src/AHv2/amplicon_hunter"
            ;;
        "AHv3")
            COMPRESS_BINARY="./src/AHv3/amplicon_hunter"
            ;;
        "AHv4")
            COMPRESS_BINARY="./src/AHv4/amplicon_hunter"
            ;;
		"AHv5")
            COMPRESS_BINARY="./src/AHv5/amplicon_hunter"
            ;;
    esac
    
    echo "Compressing genomes with ${version} (${DATASET_NAME})..."
    /usr/bin/time -v ${COMPRESS_BINARY} compress \
        --input-dir "${genome_dir}" \
        --output "${output_dir}" \
        --threads ${thread_count} \
        --batch-size 256 \
        2>&1 | tee "${LOG_FILE}"
    
    # Extract and save compression metrics
    if [ -f "${LOG_FILE}" ]; then
        echo "Compression timing saved to ${LOG_FILE}"
        
        # Extract key metrics to separate files for easy parsing
        grep "User time" "${LOG_FILE}" | awk '{print $4}' > "${LOG_FILE%.log}.user_time"
        grep "System time" "${LOG_FILE}" | awk '{print $4}' > "${LOG_FILE%.log}.sys_time"
        
        # Extract and convert elapsed time
        local elapsed_raw=$(grep "Elapsed (wall clock)" "${LOG_FILE}" | awk '{print $8}')
        if [ -n "$elapsed_raw" ]; then
            convert_time_to_seconds "$elapsed_raw" > "${LOG_FILE%.log}.real_time"
        fi
        
        grep "Maximum resident set" "${LOG_FILE}" | awk '{print $6}' > "${LOG_FILE%.log}.max_memory"
        
        # Print summary
        echo "  Compression completed in: ${elapsed_raw}"
        echo "  Max memory used: $(grep "Maximum resident set" "${LOG_FILE}" | awk '{print $6}') KB"
    fi
    
    # Clear cache after compression so it doesn't affect benchmarks
    clear_cache
}

# Function to run a single benchmark
run_benchmark() {
    local version="$1"
    local input_file="$2"
    local primers="$3"
    local output_prefix="$4"
    local mismatches="$5"
    local thread_count="${6:-$THREADS}"  # Optional thread count
    local clear_cache_flag="${7:-true}"  # Optional cache clearing flag
    
    echo "      Running ${version}..."
    
    # Add cooldown before run
    sleep ${COOLDOWN_SECONDS}
    
    # Clear cache if requested
    if [ "${clear_cache_flag}" = true ]; then
        clear_cache
    else
        echo "      Running with warm cache..."
    fi
    
    # Log system state
    echo "      CPU temps (if available):"
    sensors 2>/dev/null | grep -E "Core|Package" | head -5 || echo "      (thermal data unavailable)"
    
    local OUTPUT="${output_prefix}"
    
    case "${version}" in
        "AHv1")
            mkdir -p "${OUTPUT}_output"
            timeout ${TIMEOUT_SECONDS} /usr/bin/time -v /home/rye/miniconda3/bin/ampliconhunter run \
                "${input_file}" \
                "${primers}" \
                "${OUTPUT}_output" \
                --threads ${thread_count} \
                --mismatches ${mismatches} \
                --clamp ${CLAMP} \
                --Lmin ${MIN_LENGTH} \
                --Lmax ${MAX_LENGTH} \
                --fb-len ${FB_LEN} \
                --rb-len ${RB_LEN} \
                --trim-primers \
                --include-offtarget \
                --clobber \
                2>&1 | tee "${OUTPUT}.log"
            ;;
        "AHv2"|"AHv3"|"AHv4"|"AHv5")
            local BINARY
            case "${version}" in
                "AHv2") BINARY="./src/AHv2/amplicon_hunter" ;;
                "AHv3") BINARY="./src/AHv3/amplicon_hunter" ;;
                "AHv4") BINARY="./src/AHv4/amplicon_hunter" ;;
				"AHv5") BINARY="./src/AHv5/amplicon_hunter" ;;
            esac
            
            timeout ${TIMEOUT_SECONDS} /usr/bin/time -v ${BINARY} run \
                --input "${input_file}" \
                --primers "${primers}" \
                --output "${OUTPUT}_amplicons.fa" \
                --threads ${thread_count} \
                --mismatches ${mismatches} \
                --clamp ${CLAMP} \
                --min-length ${MIN_LENGTH} \
                --max-length ${MAX_LENGTH} \
                --fb-len ${FB_LEN} \
                --rb-len ${RB_LEN} \
                --trim-primers \
                --include-offtarget \
                2>&1 | tee "${OUTPUT}.log"
            ;;
    esac
    
    # Extract metrics with proper time conversion
    if [ $? -eq 0 ]; then
        grep "User time" "${OUTPUT}.log" | awk '{print $4}' > "${OUTPUT}.user_time"
        grep "System time" "${OUTPUT}.log" | awk '{print $4}' > "${OUTPUT}.sys_time"
        
        # Extract and convert elapsed time
        local elapsed_raw=$(grep "Elapsed (wall clock)" "${OUTPUT}.log" | awk '{print $8}')
        if [ -n "$elapsed_raw" ]; then
            convert_time_to_seconds "$elapsed_raw" > "${OUTPUT}.real_time"
        fi
        
        grep "Maximum resident set" "${OUTPUT}.log" | awk '{print $6}' > "${OUTPUT}.max_memory"
        echo "SUCCESS" > "${OUTPUT}.status"
    else
        echo "TIMEOUT_OR_ERROR" > "${OUTPUT}.status"
    fi
    
    # Kill any lingering processes
    pkill -f amplicon_hunter 2>/dev/null || true
    
    # Log completion
    echo "      Completed ${version}"
}

# Initialize system
echo "Initializing system for benchmarking..."
set_performance_governor

# TEST 1: Input Size Scaling
if [ "$RUN_TEST1" = true ]; then
    echo "=========================================="
    echo "TEST 1: INPUT SIZE SCALING"
    echo "=========================================="

    GENOME_SIZES=(G006400 G012800 G025600 G051200 G102400 G204800)
    GENOME_COUNTS=(6400 12800 25600 51200 102400 204800)

    for i in "${!GENOME_SIZES[@]}"; do
        GENOME_SIZE="${GENOME_SIZES[$i]}"
        GENOME_COUNT="${GENOME_COUNTS[$i]}"
        GENOME_DIR="${FASTA_DIR}/${GENOME_SIZE}"
        
        echo "Processing ${GENOME_COUNT} genomes..."
        
        # Prepare file list for AHv1
        find "${GENOME_DIR}" -name "*.fa" -type f > "${RESULTS_DIR}/test1/${GENOME_SIZE}_filelist.txt"
        
        # Prepare compression for versions that need it
        for version in "${VERSIONS[@]}"; do
            if [[ "${version}" == "AHv1" ]]; then
                continue  # AHv1 doesn't use compression
            fi
            if [ "${version}" == "AHv3" ]; then
                # AHv3 uses AHv4's compression
                continue
            fi
			if [ "${version}" == "AHv5" ]; then
                # AHv5 uses AHv4's compression
                continue
            fi
            compress_genomes "${version}" "${GENOME_DIR}" "${COMPRESSED_DIR}/${version}/${GENOME_SIZE}"
        done
        
        # Create list files for compressed versions
        for version in "${VERSIONS[@]}"; do
            COMPRESS_VERSION="${version}"
            if [ "${version}" == "AHv3" ]; then
                COMPRESS_VERSION="AHv4"  # AHv3 uses AHv4's compressed files
            fi
			if [ "${version}" == "AHv5" ]; then
                COMPRESS_VERSION="AHv4"  # AHv5 uses AHv4's compressed files
            fi
            
            if [ -d "${COMPRESSED_DIR}/${COMPRESS_VERSION}/${GENOME_SIZE}" ]; then
                find "$(realpath ${COMPRESSED_DIR}/${COMPRESS_VERSION}/${GENOME_SIZE})" -name "*.2bit" > "${RESULTS_DIR}/test1/${GENOME_SIZE}_${version}.list"
            fi
        done
        
        # Run benchmarks with repeats and rotation
        for repeat in $(seq 1 ${REPEATS}); do
            echo "  Repeat ${repeat}/${REPEATS} for ${GENOME_COUNT} genomes"
            
            # Rotate through all versions dynamically
            num_versions=${#VERSIONS[@]}
            ORDER=()
            for ((i=0; i<num_versions; i++)); do
                idx=$(( (i + repeat - 1) % num_versions ))
                ORDER+=("${VERSIONS[$idx]}")
            done
            
            echo "    Order for this repeat: ${ORDER[*]}"
            
            for version in "${ORDER[@]}"; do
                case $version in
                    "AHv1")
                        run_benchmark "AHv1" \
                            "${RESULTS_DIR}/test1/${GENOME_SIZE}_filelist.txt" \
                            "${PRIMERS}" \
                            "${RESULTS_DIR}/test1/AHv1/${GENOME_SIZE}_repeat${repeat}" \
                            "2"
                        ;;
                    "AHv2"|"AHv3"|"AHv4"|"AHv5")
                        LIST_FILE="${RESULTS_DIR}/test1/${GENOME_SIZE}_${version}.list"
                        if [ -f "${LIST_FILE}" ]; then
                            run_benchmark "${version}" \
                                "${LIST_FILE}" \
                                "${PRIMERS}" \
                                "${RESULTS_DIR}/test1/${version}/${GENOME_SIZE}_repeat${repeat}" \
                                "2"
                        fi
                        ;;
                esac
            done
        done
    done
fi  # End of TEST 1

# TEST 2: Primer N-base Experiment
if [ "$RUN_TEST2" = true ]; then
    echo "=========================================="
    echo "TEST 2: PRIMER N-BASE EXPERIMENT"
    echo "=========================================="

    GENOME_SIZE="G012800"
    GENOME_DIR="${FASTA_DIR}/${GENOME_SIZE}"
    N_COUNTS=(0 2 4 6)

    # Prepare file list for AHv1
    find "${GENOME_DIR}" -name "*.fa" -type f > "${RESULTS_DIR}/test2/${GENOME_SIZE}_filelist.txt"

    # Prepare compression for versions that need it
    for version in "${VERSIONS[@]}"; do
        if [[ "${version}" == "AHv1" ]]; then
            continue  # AHv1 doesn't use compression
        fi
        if [ "${version}" == "AHv3" ]; then
            continue
        fi
		if [ "${version}" == "AHv5" ]; then
            continue
        fi
        compress_genomes "${version}" "${GENOME_DIR}" "${COMPRESSED_DIR}/${version}/${GENOME_SIZE}"
    done

    # Create list files for compressed versions
    for version in "${VERSIONS[@]}"; do
        if [[ "${version}" == "AHv1" ]]; then
            continue  # AHv1 doesn't use compression
        fi
        COMPRESS_VERSION="${version}"
        if [ "${version}" == "AHv3" ]; then
            COMPRESS_VERSION="AHv4"
        fi
		if [ "${version}" == "AHv5" ]; then
            COMPRESS_VERSION="AHv4"
        fi
        
        if [ -d "${COMPRESSED_DIR}/${COMPRESS_VERSION}/${GENOME_SIZE}" ]; then
            find "$(realpath ${COMPRESSED_DIR}/${COMPRESS_VERSION}/${GENOME_SIZE})" -name "*.2bit" > "${RESULTS_DIR}/test2/${GENOME_SIZE}_${version}.list"
        fi
    done

    for n_count in "${N_COUNTS[@]}"; do
        echo "Testing with ${n_count} N bases at 3' end..."
        
        # Create modified primers file
        PRIMER_FILE="${RESULTS_DIR}/test2/primers_${n_count}N.txt"
        if [ ${n_count} -eq 0 ]; then
            cp "${PRIMERS}" "${PRIMER_FILE}"
        else
            # Read original primers
            FWD_PRIMER=$(head -1 "${PRIMERS}" | cut -f2)
            REV_PRIMER=$(tail -1 "${PRIMERS}" | cut -f2)
            
            # Add N bases to 3' end
            N_STRING=$(printf 'N%.0s' $(seq 1 ${n_count}))
            FWD_MODIFIED="${FWD_PRIMER}${N_STRING}"
            REV_MODIFIED="${REV_PRIMER}${N_STRING}"
            
            # Write modified primers
            echo -e "${FWD_MODIFIED}\n${REV_MODIFIED}" > "${PRIMER_FILE}"
        fi
        
        # Run benchmarks with repeats and rotation
        for repeat in $(seq 1 ${REPEATS}); do
            echo "  Repeat ${repeat}/${REPEATS} for ${n_count}N"
            
            # Rotate through all versions dynamically
            num_versions=${#VERSIONS[@]}
            ORDER=()
            for ((i=0; i<num_versions; i++)); do
                idx=$(( (i + repeat - 1) % num_versions ))
                ORDER+=("${VERSIONS[$idx]}")
            done
            
            echo "    Order for this repeat: ${ORDER[*]}"
            
            for version in "${ORDER[@]}"; do
                case $version in
                    "AHv1")
                        run_benchmark "AHv1" \
                            "${RESULTS_DIR}/test2/${GENOME_SIZE}_filelist.txt" \
                            "${PRIMER_FILE}" \
                            "${RESULTS_DIR}/test2/AHv1/${n_count}N_repeat${repeat}" \
                            "2"
                        ;;
                    "AHv2"|"AHv3"|"AHv4"|"AHv5")
                        LIST_FILE="${RESULTS_DIR}/test2/${GENOME_SIZE}_${version}.list"
                        if [ -f "${LIST_FILE}" ]; then
                            run_benchmark "${version}" \
                                "${LIST_FILE}" \
                                "${PRIMER_FILE}" \
                                "${RESULTS_DIR}/test2/${version}/${n_count}N_repeat${repeat}" \
                                "2"
                        fi
                        ;;
                esac
            done
        done
    done
fi  # End of TEST 2

# TEST 3: Mismatch Variation
if [ "$RUN_TEST3" = true ]; then
    echo "=========================================="
    echo "TEST 3: MISMATCH VARIATION"
    echo "=========================================="

    GENOME_SIZE="G012800"
    MISMATCH_COUNTS=(0 1 2 3 4 5 6)

    for mismatch in "${MISMATCH_COUNTS[@]}"; do
        echo "Testing with ${mismatch} mismatches..."
        
        # Run benchmarks with repeats and rotation
        for repeat in $(seq 1 ${REPEATS}); do
            echo "  Repeat ${repeat}/${REPEATS} for ${mismatch} mismatches"
            
            # Rotate through all versions dynamically
            num_versions=${#VERSIONS[@]}
            ORDER=()
            for ((i=0; i<num_versions; i++)); do
                idx=$(( (i + repeat - 1) % num_versions ))
                ORDER+=("${VERSIONS[$idx]}")
            done
            
            echo "    Order for this repeat: ${ORDER[*]}"
            
            for version in "${ORDER[@]}"; do
                case $version in
                    "AHv1")
                        run_benchmark "AHv1" \
                            "${RESULTS_DIR}/test2/${GENOME_SIZE}_filelist.txt" \
                            "${PRIMERS}" \
                            "${RESULTS_DIR}/test3/AHv1/mm${mismatch}_repeat${repeat}" \
                            "${mismatch}"
                        ;;
                    "AHv2"|"AHv3"|"AHv4"|"AHv5")
                        LIST_FILE="${RESULTS_DIR}/test2/${GENOME_SIZE}_${version}.list"
                        if [ -f "${LIST_FILE}" ]; then
                            run_benchmark "${version}" \
                                "${LIST_FILE}" \
                                "${PRIMERS}" \
                                "${RESULTS_DIR}/test3/${version}/mm${mismatch}_repeat${repeat}" \
                                "${mismatch}"
                        fi
                        ;;
                esac
            done
        done
    done
fi  # End of TEST 3

# TEST 4: Thread scaling
if [ "$RUN_TEST4" = true ]; then
    echo "=========================================="
    echo "TEST 4: THREAD SCALING"
    echo "=========================================="

    # Use G012800 dataset for thread scaling test
    THREAD_DATASET="G012800"
    THREAD_GENOME_DIR="${FASTA_DIR}/${THREAD_DATASET}"
    THREAD_GENOME_COUNT=12800

    # Thread counts to test
    THREAD_COUNTS=(1 2 4 8 16 32 64 96 128 160 190)

    # Create file list for AHv1
    THREAD_FILELIST="${RESULTS_DIR}/test4/${THREAD_DATASET}_filelist.txt"
    if [ ! -f "${THREAD_FILELIST}" ]; then
        echo "Creating file list for thread scaling test..."
        find "${THREAD_GENOME_DIR}" -name "*.fa" -type f | head -${THREAD_GENOME_COUNT} | sort > "${THREAD_FILELIST}"
    fi

    # Test each thread count
    for thread_count in "${THREAD_COUNTS[@]}"; do
        echo ""
        echo "Testing with ${thread_count} threads..."
        
        # For compressed versions, recompress with the specified thread count
        for version in "${VERSIONS[@]}"; do
            if [ "${version}" == "AHv3" ]; then
                continue  # AHv3 uses AHv4's compression
            fi
			if [ "${version}" == "AHv5" ]; then
                continue  # AHv5 uses AHv4's compression
            fi
            
            echo "  Compressing ${THREAD_GENOME_COUNT} genomes for ${version} with ${thread_count} threads..."
            COMPRESS_DIR="${COMPRESSED_DIR}/${version}/${THREAD_DATASET}_t${thread_count}"
            
            if [ ! -d "${COMPRESS_DIR}" ] || [ -z "$(ls -A ${COMPRESS_DIR} 2>/dev/null)" ]; then
                compress_genomes "${version}" "${THREAD_GENOME_DIR}" "${COMPRESS_DIR}" "${thread_count}"
            else
                echo "    Using existing compressed data in ${COMPRESS_DIR}"
            fi
        done
        
        # Create list files for compressed versions
        for version in "${VERSIONS[@]}"; do
            COMPRESS_VERSION="${version}"
            if [ "${version}" == "AHv3" ]; then
                COMPRESS_VERSION="AHv4"
            fi
			if [ "${version}" == "AHv5" ]; then
                COMPRESS_VERSION="AHv4"
            fi
            
            COMPRESS_DIR="${COMPRESSED_DIR}/${COMPRESS_VERSION}/${THREAD_DATASET}_t${thread_count}"
            THREAD_LIST="${RESULTS_DIR}/test4/${THREAD_DATASET}_t${thread_count}_${version}.list"
            find "$(realpath ${COMPRESS_DIR})" -name "*.2bit" | sort > "${THREAD_LIST}"
        done
        
        # Run 3 repeats for each thread count
        for repeat in $(seq 1 ${REPEATS}); do
            echo "  Repeat ${repeat}/${REPEATS} for ${thread_count} threads"
            
            # Rotate through all versions dynamically
            num_versions=${#VERSIONS[@]}
            ORDER=()
            for ((i=0; i<num_versions; i++)); do
                idx=$(( (i + repeat - 1) % num_versions ))
                ORDER+=("${VERSIONS[$idx]}")
            done
            
            echo "    Order for this repeat: ${ORDER[*]}"
            
            for version in "${ORDER[@]}"; do
                echo "      Running ${version} with ${thread_count} threads..."
                
                case $version in
                    "AHv1")
                        run_benchmark "AHv1" \
                            "${THREAD_FILELIST}" \
                            "${PRIMERS}" \
                            "${RESULTS_DIR}/test4/AHv1/t${thread_count}_repeat${repeat}" \
                            "2" \
                            "${thread_count}"
                        ;;
                    "AHv2"|"AHv3"|"AHv4"|"AHv5")
                        THREAD_LIST="${RESULTS_DIR}/test4/${THREAD_DATASET}_t${thread_count}_${version}.list"
                        if [ -f "${THREAD_LIST}" ]; then
                            run_benchmark "${version}" \
                                "${THREAD_LIST}" \
                                "${PRIMERS}" \
                                "${RESULTS_DIR}/test4/${version}/t${thread_count}_repeat${repeat}" \
                                "2" \
                                "${thread_count}"
                        fi
                        ;;
                esac
            done
            
            # Cooldown between runs
            sleep ${COOLDOWN_SECONDS}
        done
    done
fi  # End of TEST 4

# TEST 5: Hot vs Cold Cache
if [ "$RUN_TEST5" = true ]; then
    echo "=========================================="
    echo "TEST 5: HOT VS COLD CACHE"
    echo "=========================================="

    CACHE_DATASET="G012800"
    CACHE_GENOME_DIR="${FASTA_DIR}/${CACHE_DATASET}"
    
    # Create file lists (reuse from test2 if available)
    CACHE_FILELIST="${RESULTS_DIR}/test5/${CACHE_DATASET}_filelist.txt"
    if [ ! -f "${CACHE_FILELIST}" ]; then
        find "${CACHE_GENOME_DIR}" -name "*.fa" -type f > "${CACHE_FILELIST}"
    fi
    
    # Create list files for compressed versions (reuse from test2)
    for version in "${VERSIONS[@]}"; do
        if [[ "${version}" == "AHv1" ]]; then
            continue  # AHv1 doesn't use compression
        fi
        SOURCE_LIST="${RESULTS_DIR}/test2/${CACHE_DATASET}_${version}.list"
        DEST_LIST="${RESULTS_DIR}/test5/${CACHE_DATASET}_${version}.list"
        if [ -f "${SOURCE_LIST}" ]; then
            cp "${SOURCE_LIST}" "${DEST_LIST}"
        fi
    done
    
    echo "Running cold cache tests (with cache clearing)..."
    for version in "${VERSIONS[@]}"; do
        echo "  Testing ${version} with cold cache..."
        
        for repeat in $(seq 1 ${REPEATS}); do
            echo "    Cold cache repeat ${repeat}/${REPEATS}"
            
            case $version in
                "AHv1")
                    run_benchmark "AHv1" \
                        "${CACHE_FILELIST}" \
                        "${PRIMERS}" \
                        "${RESULTS_DIR}/test5/AHv1/cold_repeat${repeat}" \
                        "2" \
                        "${THREADS}" \
                        true  # Clear cache
                    ;;
                "AHv2"|"AHv3"|"AHv4"|"AHv5")
                    LIST_FILE="${RESULTS_DIR}/test5/${CACHE_DATASET}_${version}.list"
                    if [ -f "${LIST_FILE}" ]; then
                        run_benchmark "${version}" \
                            "${LIST_FILE}" \
                            "${PRIMERS}" \
                            "${RESULTS_DIR}/test5/${version}/cold_repeat${repeat}" \
                            "2" \
                            "${THREADS}" \
                            true  # Clear cache
                    fi
                    ;;
            esac
        done
    done
    
    echo ""
    echo "Running hot cache tests (without cache clearing)..."
    for version in "${VERSIONS[@]}"; do
        echo "  Testing ${version} with hot cache..."
        
        # Warm up the cache with one run
        echo "    Warming up cache for ${version}..."
        case $version in
            "AHv1")
                run_benchmark "AHv1" \
                    "${CACHE_FILELIST}" \
                    "${PRIMERS}" \
                    "${RESULTS_DIR}/test5/AHv1/warmup" \
                    "2" \
                    "${THREADS}" \
                    true  # Clear cache first for fair warmup
                ;;
            "AHv2"|"AHv3"|"AHv4"|"AHv5")
                LIST_FILE="${RESULTS_DIR}/test5/${CACHE_DATASET}_${version}.list"
                if [ -f "${LIST_FILE}" ]; then
                    run_benchmark "${version}" \
                        "${LIST_FILE}" \
                        "${PRIMERS}" \
                        "${RESULTS_DIR}/test5/${version}/warmup" \
                        "2" \
                        "${THREADS}" \
                        true  # Clear cache first for fair warmup
                fi
                ;;
        esac
        
        # Now run hot cache tests (no cache clearing)
        for repeat in $(seq 1 ${REPEATS}); do
            echo "    Hot cache repeat ${repeat}/${REPEATS}"
            
            case $version in
                "AHv1")
                    run_benchmark "AHv1" \
                        "${CACHE_FILELIST}" \
                        "${PRIMERS}" \
                        "${RESULTS_DIR}/test5/AHv1/hot_repeat${repeat}" \
                        "2" \
                        "${THREADS}" \
                        false  # Don't clear cache
                    ;;
                "AHv2"|"AHv3"|"AHv4"|"AHv5")
                    LIST_FILE="${RESULTS_DIR}/test5/${CACHE_DATASET}_${version}.list"
                    if [ -f "${LIST_FILE}" ]; then
                        run_benchmark "${version}" \
                            "${LIST_FILE}" \
                            "${PRIMERS}" \
                            "${RESULTS_DIR}/test5/${version}/hot_repeat${repeat}" \
                            "2" \
                            "${THREADS}" \
                            false  # Don't clear cache
                    fi
                    ;;
            esac
        done
    done
fi  # End of TEST 5

# TEST 6: Primer Pair Comparison
if [ "$RUN_TEST6" = true ]; then
    echo "=========================================="
    echo "TEST 6: PRIMER PAIR COMPARISON"
    echo "=========================================="

    GENOME_SIZE="G012800"
    GENOME_DIR="${FASTA_DIR}/${GENOME_SIZE}"
    
    # Three primer pairs to test
    PRIMER_PAIRS=("V3V4" "Titan" "V1V9")
    PRIMER_FILES=(
        "${DATA_DIR}/primers/primers_V3V4.txt"
        "${DATA_DIR}/primers/primers_Titan.txt"
        "${DATA_DIR}/primers/primers_V1V9.txt"
    )
    
    # Prepare file list for AHv1 (reuse from test2 if available)
    if [ -f "${RESULTS_DIR}/test2/${GENOME_SIZE}_filelist.txt" ]; then
        cp "${RESULTS_DIR}/test2/${GENOME_SIZE}_filelist.txt" "${RESULTS_DIR}/test6/${GENOME_SIZE}_filelist.txt"
    else
        find "${GENOME_DIR}" -name "*.fa" -type f > "${RESULTS_DIR}/test6/${GENOME_SIZE}_filelist.txt"
    fi
    
    # Create list files for compressed versions (reuse from test2)
    for version in "${VERSIONS[@]}"; do
        if [[ "${version}" == "AHv1" ]]; then
            continue  # AHv1 doesn't use compression
        fi
        SOURCE_LIST="${RESULTS_DIR}/test2/${GENOME_SIZE}_${version}.list"
        DEST_LIST="${RESULTS_DIR}/test6/${GENOME_SIZE}_${version}.list"
        if [ -f "${SOURCE_LIST}" ]; then
            cp "${SOURCE_LIST}" "${DEST_LIST}"
        else
            # If test2 wasn't run, create the list files
            COMPRESS_VERSION="${version}"
            if [ "${version}" == "AHv3" ]; then
                COMPRESS_VERSION="AHv4"  # AHv3 uses AHv4's compressed files
            fi
			if [ "${version}" == "AHv5" ]; then
                COMPRESS_VERSION="AHv4"  # AHv5 uses AHv4's compressed files
            fi
            
            if [ -d "${COMPRESSED_DIR}/${COMPRESS_VERSION}/${GENOME_SIZE}" ]; then
                find "$(realpath ${COMPRESSED_DIR}/${COMPRESS_VERSION}/${GENOME_SIZE})" -name "*.2bit" > "${DEST_LIST}"
            fi
        fi
    done
    
    # Test each primer pair
    for i in "${!PRIMER_PAIRS[@]}"; do
        PRIMER_NAME="${PRIMER_PAIRS[$i]}"
        PRIMER_FILE="${PRIMER_FILES[$i]}"
        
        if [ ! -f "${PRIMER_FILE}" ]; then
            echo "Warning: Primer file ${PRIMER_FILE} not found, skipping ${PRIMER_NAME}"
            continue
        fi
        
        echo "Testing with ${PRIMER_NAME} primers..."
        
        # Run benchmarks with repeats and rotation
        for repeat in $(seq 1 ${REPEATS}); do
            echo "  Repeat ${repeat}/${REPEATS} for ${PRIMER_NAME} primers"
            
            # Rotate through all versions dynamically
            num_versions=${#VERSIONS[@]}
            ORDER=()
            for ((i=0; i<num_versions; i++)); do
                idx=$(( (i + repeat - 1) % num_versions ))
                ORDER+=("${VERSIONS[$idx]}")
            done
            
            echo "    Order for this repeat: ${ORDER[*]}"
            
            for version in "${ORDER[@]}"; do
                case $version in
                    "AHv1")
                        run_benchmark "AHv1" \
                            "${RESULTS_DIR}/test6/${GENOME_SIZE}_filelist.txt" \
                            "${PRIMER_FILE}" \
                            "${RESULTS_DIR}/test6/AHv1/${PRIMER_NAME}_repeat${repeat}" \
                            "2"
                        ;;
                    "AHv2"|"AHv3"|"AHv4"|"AHv5")
                        LIST_FILE="${RESULTS_DIR}/test6/${GENOME_SIZE}_${version}.list"
                        if [ -f "${LIST_FILE}" ]; then
                            run_benchmark "${version}" \
                                "${LIST_FILE}" \
                                "${PRIMER_FILE}" \
                                "${RESULTS_DIR}/test6/${version}/${PRIMER_NAME}_repeat${repeat}" \
                                "2"
                        fi
                        ;;
                esac
            done
        done
    done
fi  # End of TEST 6

echo "Benchmarking complete! Generating report..."

# Generate Python report script for all versions and all tests
cat > "${RESULTS_DIR}/generate_report.py" << 'PYEOF'
import os
import sys
import json
import glob
import statistics

results_dir = sys.argv[1] if len(sys.argv) > 1 else "results"

def read_metric(filepath):
    """Read a single metric value from file"""
    try:
        with open(filepath, 'r') as f:
            content = f.read().strip()
            if content == "TIMEOUT_OR_ERROR":
                return None
            # Handle time format if present
            if ':' in content:
                parts = content.split(':')
                if len(parts) == 3:  # H:MM:SS
                    return float(parts[0]) * 3600 + float(parts[1]) * 60 + float(parts[2])
                elif len(parts) == 2:  # M:SS
                    return float(parts[0]) * 60 + float(parts[1])
            return float(content) if content else None
    except:
        return None

def collect_test1_results(results_dir):
    """Collect results for Test 1: Input Size Scaling"""
    results = {}
    genome_sizes = ['G006400', 'G012800', 'G025600', 'G051200', 'G102400', 'G204800']
    genome_counts = [6400, 12800, 25600, 51200, 102400, 204800]
    
    for version in ['AHv1', 'AHv2', 'AHv3', 'AHv4', 'AHv5']:
        results[version] = {}
        for size, count in zip(genome_sizes, genome_counts):
            metrics = {'real_time': [], 'user_time': [], 'sys_time': [], 'max_memory': []}
            
            for repeat in range(1, 6):
                base_path = f"{results_dir}/test1/{version}/{size}_repeat{repeat}"
                
                real_time = read_metric(f"{base_path}.real_time")
                user_time = read_metric(f"{base_path}.user_time")
                sys_time = read_metric(f"{base_path}.sys_time")
                max_memory = read_metric(f"{base_path}.max_memory")
                
                if real_time is not None:
                    metrics['real_time'].append(real_time)
                if user_time is not None:
                    metrics['user_time'].append(user_time)
                if sys_time is not None:
                    metrics['sys_time'].append(sys_time)
                if max_memory is not None:
                    metrics['max_memory'].append(max_memory)
            
            results[version][count] = metrics
    
    return results

def collect_test2_results(results_dir):
    """Collect results for Test 2: Primer N-base Experiment"""
    results = {}
    n_counts = [0, 2, 4, 6]
    
    for version in ['AHv1', 'AHv2', 'AHv3', 'AHv4', 'AHv5']:
        results[version] = {}
        for n_count in n_counts:
            metrics = {'real_time': [], 'user_time': [], 'sys_time': [], 'max_memory': []}
            
            for repeat in range(1, 6):
                base_path = f"{results_dir}/test2/{version}/{n_count}N_repeat{repeat}"
                
                real_time = read_metric(f"{base_path}.real_time")
                user_time = read_metric(f"{base_path}.user_time")
                sys_time = read_metric(f"{base_path}.sys_time")
                max_memory = read_metric(f"{base_path}.max_memory")
                
                if real_time is not None:
                    metrics['real_time'].append(real_time)
                if user_time is not None:
                    metrics['user_time'].append(user_time)
                if sys_time is not None:
                    metrics['sys_time'].append(sys_time)
                if max_memory is not None:
                    metrics['max_memory'].append(max_memory)
            
            results[version][f"{n_count}N"] = metrics
    
    return results

def collect_test3_results(results_dir):
    """Collect results for Test 3: Mismatch Variation"""
    results = {}
    mismatch_counts = [0, 1, 2, 3, 4, 5, 6]
    
    for version in ['AHv1', 'AHv2', 'AHv3', 'AHv4', 'AHv5']:
        results[version] = {}
        for mm in mismatch_counts:
            metrics = {'real_time': [], 'user_time': [], 'sys_time': [], 'max_memory': []}
            
            for repeat in range(1, 6):
                base_path = f"{results_dir}/test3/{version}/mm{mm}_repeat{repeat}"
                
                real_time = read_metric(f"{base_path}.real_time")
                user_time = read_metric(f"{base_path}.user_time")
                sys_time = read_metric(f"{base_path}.sys_time")
                max_memory = read_metric(f"{base_path}.max_memory")
                
                if real_time is not None:
                    metrics['real_time'].append(real_time)
                if user_time is not None:
                    metrics['user_time'].append(user_time)
                if sys_time is not None:
                    metrics['sys_time'].append(sys_time)
                if max_memory is not None:
                    metrics['max_memory'].append(max_memory)
            
            results[version][mm] = metrics
    
    return results

def collect_test4_results(results_dir):
    """Collect results for Test 4: Thread Scaling"""
    results = {}
    thread_counts = [1, 2, 4, 8, 16, 32, 64, 96, 128, 160, 190]
    
    for version in ['AHv1', 'AHv2', 'AHv3', 'AHv4', 'AHv5']:
        results[version] = {}
        for threads in thread_counts:
            metrics = {'real_time': [], 'user_time': [], 'sys_time': [], 'max_memory': []}
            
            for repeat in range(1, 6):
                base_path = f"{results_dir}/test4/{version}/t{threads}_repeat{repeat}"
                
                real_time = read_metric(f"{base_path}.real_time")
                user_time = read_metric(f"{base_path}.user_time")
                sys_time = read_metric(f"{base_path}.sys_time")
                max_memory = read_metric(f"{base_path}.max_memory")
                
                if real_time is not None:
                    metrics['real_time'].append(real_time)
                if user_time is not None:
                    metrics['user_time'].append(user_time)
                if sys_time is not None:
                    metrics['sys_time'].append(sys_time)
                if max_memory is not None:
                    metrics['max_memory'].append(max_memory)
            
            results[version][threads] = metrics
    
    return results

def collect_test5_results(results_dir):
    """Collect results for Test 5: Hot vs Cold Cache"""
    results = {}
    
    for version in ['AHv1', 'AHv2', 'AHv3', 'AHv4', 'AHv5']:
        results[version] = {'cold': {'real_time': [], 'user_time': [], 'sys_time': [], 'max_memory': []},
                           'hot': {'real_time': [], 'user_time': [], 'sys_time': [], 'max_memory': []}}
        
        # Cold cache results
        for repeat in range(1, 6):
            base_path = f"{results_dir}/test5/{version}/cold_repeat{repeat}"
            
            real_time = read_metric(f"{base_path}.real_time")
            user_time = read_metric(f"{base_path}.user_time")
            sys_time = read_metric(f"{base_path}.sys_time")
            max_memory = read_metric(f"{base_path}.max_memory")
            
            if real_time is not None:
                results[version]['cold']['real_time'].append(real_time)
            if user_time is not None:
                results[version]['cold']['user_time'].append(user_time)
            if sys_time is not None:
                results[version]['cold']['sys_time'].append(sys_time)
            if max_memory is not None:
                results[version]['cold']['max_memory'].append(max_memory)
        
        # Hot cache results
        for repeat in range(1, 6):
            base_path = f"{results_dir}/test5/{version}/hot_repeat{repeat}"
            
            real_time = read_metric(f"{base_path}.real_time")
            user_time = read_metric(f"{base_path}.user_time")
            sys_time = read_metric(f"{base_path}.sys_time")
            max_memory = read_metric(f"{base_path}.max_memory")
            
            if real_time is not None:
                results[version]['hot']['real_time'].append(real_time)
            if user_time is not None:
                results[version]['hot']['user_time'].append(user_time)
            if sys_time is not None:
                results[version]['hot']['sys_time'].append(sys_time)
            if max_memory is not None:
                results[version]['hot']['max_memory'].append(max_memory)
    
    return results

def collect_test6_results(results_dir):
    """Collect results for Test 6: Primer Pair Comparison"""
    results = {}
    primer_pairs = ['V3V4', 'Titan', 'V1V9']
    
    for version in ['AHv1', 'AHv2', 'AHv3', 'AHv4', 'AHv5']:
        results[version] = {}
        for primer in primer_pairs:
            metrics = {'real_time': [], 'user_time': [], 'sys_time': [], 'max_memory': []}
            
            for repeat in range(1, 6):
                base_path = f"{results_dir}/test6/{version}/{primer}_repeat{repeat}"
                
                real_time = read_metric(f"{base_path}.real_time")
                user_time = read_metric(f"{base_path}.user_time")
                sys_time = read_metric(f"{base_path}.sys_time")
                max_memory = read_metric(f"{base_path}.max_memory")
                
                if real_time is not None:
                    metrics['real_time'].append(real_time)
                if user_time is not None:
                    metrics['user_time'].append(user_time)
                if sys_time is not None:
                    metrics['sys_time'].append(sys_time)
                if max_memory is not None:
                    metrics['max_memory'].append(max_memory)
            
            results[version][primer] = metrics
    
    return results

def calculate_stats(values):
    """Calculate mean and standard deviation"""
    if not values:
        return None, None
    if len(values) == 1:
        return values[0], 0
    return statistics.mean(values), statistics.stdev(values)

def format_time(seconds):
    """Format time in seconds to human readable format"""
    if seconds is None:
        return "N/A"
    if seconds >= 3600:
        return f"{seconds/3600:.2f}h"
    elif seconds >= 60:
        return f"{seconds/60:.2f}m"
    else:
        return f"{seconds:.2f}s"

def format_memory(kb):
    """Format memory in KB to human readable format"""
    if kb is None:
        return "N/A"
    if kb >= 1048576:
        return f"{kb/1048576:.2f}GB"
    elif kb >= 1024:
        return f"{kb/1024:.2f}MB"
    else:
        return f"{kb:.2f}KB"

# Generate report
print("=" * 80)
print("AMPLICONHUNTER COMPREHENSIVE BENCHMARKING REPORT")
print("=" * 80)
print()

# Test 1 Results
print("TEST 1: INPUT SIZE SCALING")
print("-" * 40)
test1_results = collect_test1_results(results_dir)

print("\nReal Time (wall clock):")
print("Genomes\t\tAHv1\t\tAHv2\t\tAHv3\t\tAHv4\tAHv5")
for count in [6400, 12800, 25600, 51200, 102400, 204800]:
    row = f"{count}\t\t"
    for version in ['AHv1', 'AHv2', 'AHv3', 'AHv4', 'AHv5']:
        if count in test1_results[version]:
            mean, std = calculate_stats(test1_results[version][count]['real_time'])
            if mean is not None:
                row += f"{format_time(mean)}\t"
            else:
                row += "TIMEOUT\t"
        else:
            row += "N/A\t"
    print(row)

print("\nMaximum Memory Usage:")
print("Genomes\t\tAHv1\t\tAHv2\t\tAHv3\t\tAHv4\tAHv5")
for count in [6400, 12800, 25600, 51200, 102400, 204800]:
    row = f"{count}\t\t"
    for version in ['AHv1', 'AHv2', 'AHv3', 'AHv4', 'AHv5']:
        if count in test1_results[version]:
            mean, std = calculate_stats(test1_results[version][count]['max_memory'])
            if mean is not None:
                row += f"{format_memory(mean)}\t"
            else:
                row += "N/A\t"
        else:
            row += "N/A\t"
    print(row)

# Test 2 Results
print("\n" + "=" * 80)
print("TEST 2: PRIMER N-BASE EXPERIMENT")
print("-" * 40)
test2_results = collect_test2_results(results_dir)

print("\nReal Time (wall clock):")
print("N-bases\t\tAHv1\t\tAHv2\t\tAHv3\t\tAHv4\tAHv5")
for n in ['0N', '2N', '4N', '6N']:
    row = f"{n}\t\t"
    for version in ['AHv1', 'AHv2', 'AHv3', 'AHv4', 'AHv5']:
        if n in test2_results[version]:
            mean, std = calculate_stats(test2_results[version][n]['real_time'])
            if mean is not None:
                row += f"{format_time(mean)}\t"
            else:
                row += "TIMEOUT\t"
        else:
            row += "N/A\t"
    print(row)

# Test 3 Results
print("\n" + "=" * 80)
print("TEST 3: MISMATCH VARIATION")
print("-" * 40)
test3_results = collect_test3_results(results_dir)

print("\nReal Time (wall clock):")
print("Mismatches\tAHv1\t\tAHv2\t\tAHv3\t\tAHv4\tAHv5")
for mm in [0, 1, 2, 3, 4, 5, 6]:
    row = f"{mm}\t\t"
    for version in ['AHv1', 'AHv2', 'AHv3', 'AHv4', 'AHv5']:
        if mm in test3_results[version]:
            mean, std = calculate_stats(test3_results[version][mm]['real_time'])
            if mean is not None:
                row += f"{format_time(mean)}\t"
            else:
                row += "TIMEOUT\t"
        else:
            row += "N/A\t"
    print(row)

# Test 4 Results
print("\n" + "=" * 80)
print("TEST 4: THREAD SCALING")
print("-" * 40)
test4_results = collect_test4_results(results_dir)

print("\nReal Time for 190 threads:")
print("Version\t\tRuntime")
for version in ['AHv1', 'AHv2', 'AHv3', 'AHv4', 'AHv5']:
    if 190 in test4_results[version]:
        mean, std = calculate_stats(test4_results[version][190]['real_time'])
        if mean is not None:
            print(f"{version}\t\t{format_time(mean)} ± {format_time(std)}")
        else:
            print(f"{version}\t\tN/A")

# Test 5 Results
print("\n" + "=" * 80)
print("TEST 5: HOT VS COLD CACHE")
print("-" * 40)
test5_results = collect_test5_results(results_dir)

print("\nCache Performance Comparison:")
print("Version\t\tCold Cache\tHot Cache\tSpeedup")
for version in ['AHv1', 'AHv2', 'AHv3', 'AHv4', 'AHv5']:
    cold_mean, cold_std = calculate_stats(test5_results[version]['cold']['real_time'])
    hot_mean, hot_std = calculate_stats(test5_results[version]['hot']['real_time'])
    
    if cold_mean and hot_mean:
        speedup = cold_mean / hot_mean
        print(f"{version}\t\t{format_time(cold_mean)}\t{format_time(hot_mean)}\t{speedup:.2f}x")
    else:
        print(f"{version}\t\tN/A\t\tN/A\t\tN/A")

# Test 6 Results
print("\n" + "=" * 80)
print("TEST 6: PRIMER PAIR COMPARISON")
print("-" * 40)
test6_results = collect_test6_results(results_dir)

print("\nReal Time (wall clock):")
print("Primer\t\tAHv1\t\tAHv2\t\tAHv3\t\tAHv4\tAHv5")
for primer in ['V3V4', 'Titan', 'V1V9']:
    row = f"{primer}\t\t"
    for version in ['AHv1', 'AHv2', 'AHv3', 'AHv4', 'AHv5']:
        if primer in test6_results[version]:
            mean, std = calculate_stats(test6_results[version][primer]['real_time'])
            if mean is not None:
                row += f"{format_time(mean)}\t"
            else:
                row += "TIMEOUT\t"
        else:
            row += "N/A\t"
    print(row)

print("\nMaximum Memory Usage:")
print("Primer\t\tAHv1\t\tAHv2\t\tAHv3\t\tAHv4\tAHv5")
for primer in ['V3V4', 'Titan', 'V1V9']:
    row = f"{primer}\t\t"
    for version in ['AHv1', 'AHv2', 'AHv3', 'AHv4', 'AHv5']:
        if primer in test6_results[version]:
            mean, std = calculate_stats(test6_results[version][primer]['max_memory'])
            if mean is not None:
                row += f"{format_memory(mean)}\t"
            else:
                row += "N/A\t"
        else:
            row += "N/A\t"
    print(row)

print("\n" + "=" * 80)
print("Report generation complete!")

# Save detailed results to JSON
all_results = {
    'test1': test1_results,
    'test2': test2_results,
    'test3': test3_results,
    'test4': test4_results,
    'test5': test5_results,
    'test6': test6_results
}

with open(f"{results_dir}/benchmark_results.json", 'w') as f:
    json.dump(all_results, f, indent=2)

print(f"\nDetailed results saved to {results_dir}/benchmark_results.json")

PYEOF

python3 "${RESULTS_DIR}/generate_report.py" "${RESULTS_DIR}"

# Generate compression timing summary
echo ""
echo "=========================================="
echo "COMPRESSION TIMING SUMMARY"
echo "=========================================="

if [ -d "${RESULTS_DIR}/compression_timing" ]; then
    echo "Dataset          Version    Real Time    Memory"
    echo "-------          -------    ---------    ------"
    
    for log_file in ${RESULTS_DIR}/compression_timing/*_*.real_time; do
        if [ -f "$log_file" ]; then
            base_name=$(basename "$log_file" .real_time)
            version=$(echo "$base_name" | cut -d_ -f1)
            dataset=$(echo "$base_name" | cut -d_ -f2-)
            
            real_time=$(cat "$log_file" 2>/dev/null || echo "N/A")
            mem_file="${log_file%.real_time}.max_memory"
            memory=$(cat "$mem_file" 2>/dev/null || echo "N/A")
            
            if [ "$memory" != "N/A" ]; then
                # Convert KB to human readable
                memory_mb=$(echo "scale=2; $memory / 1024" | bc)
                memory="${memory_mb} MB"
            fi
            
            printf "%-15s  %-7s  %8ss    %s\n" "$dataset" "$version" "$real_time" "$memory"
        fi
    done
    echo ""
fi

echo "Benchmarking suite complete!"
echo ""
echo "SUMMARY:"
echo "- Testing all 5 versions: AHv1, AHv2, AHv3, AHv4, and AHv5"
echo "- AHv3-5 share compressed files"
echo "- Rotation pattern ensures fair testing order"
echo "- Using CPU core count: ${THREADS}"
echo "- Results saved to: ${RESULTS_DIR}"
echo ""
echo "To generate publication figures, run:"
echo "  python3 generate_publication_figures.py ${RESULTS_DIR}"
