#!/bin/tcsh
#$ -o "/gpfs0/tals/projects/Analysis/Oded_Project/output_assembly_run_test"
#$ -e "/gpfs0/tals/projects/Analysis/Oded_Project/error_assembly_run_test"

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <arg1> <arg2>" >> "/gpfs0/tals/projects/Analysis/Oded_Project/error_assembly_test_$1"
    exit 1
    
fi

# Assign arguments to variables
arg1="$1" # accession number
echo "$arg1"
arg2="$2" # links to fastq files separated with ';'
echo "$arg2"

# Define the directory path
dir="../Input/$arg1"

# Create the directory if it doesn't exist
if [ ! -d "$dir" ]; then
    mkdir -p "$dir"
fi

# Loop through the links in arg2 and download the files only if absent or size < expected
IFS=';' read -ra links <<< "$arg2"
for link in "${links[@]}"; do
    filename=$(basename "$link")
    output_path="$dir/$filename"
    echo "Working with $link" >> /gpfs0/tals/projects/Analysis/Oded_Project/output_assembly_test_$argv[1]

    # Check if the file is absent or its size is less than expected

    tmp_file=$(mktemp /tmp/wget_output.XXXXXX)
    wget --spider "$link" &> "$tmp_file"
    expected_size=$(grep 'Length:' "$tmp_file" | awk '{print $2}')
    echo "Expected size of file is $expected_size bytes"

    # Clean up temporary file
    rm "$tmp_file"

    # Check actual file size
    if [ ! -f "$output_path" ]; then
        actual_size=0
    else
        actual_size=$(stat -c%s "$output_path")
    fi

    echo "Actual file size is $actual_size bytes"

    # Retry the download if necessary
    if [ ! -f "$output_path" ] || [ "$actual_size" -lt "$expected_size" ]; then
        j=1
        while [ $j -le 3 ]; do
            wget -c -q "$link" -O "$output_path"

            actual_size=$(stat -c%s "$output_path")

            if [ "$actual_size" -eq "$expected_size" ]; then
                echo "Download successful: $link" >> /gpfs0/tals/projects/Analysis/Oded_Project/output_assembly_test_$arg1
                break
            elif [ $j -lt 3 ]; then
                ((j++))
                echo "Download of $link was interrupted. Retrying, attempt $j" >> /gpfs0/tals/projects/Analysis/Oded_Project/output_assembly_test_$arg1
            else
                echo "Error: Download of $link failed." >> /gpfs0/tals/projects/Analysis/Oded_Project/error_assembly_test_$arg1
                exit 1
            fi
        done
    else
        echo "Skipping download for $link as the file already exists and has the expected size." >> /gpfs0/tals/projects/Analysis/Oded_Project/output_assembly_test_$arg1
    fi

done

# Run the Python script with arguments
python /gpfs0/tals/projects/Analysis/Oded_Project/scripts/PlasmidomeAssembly.py "$arg1" "$dir"
