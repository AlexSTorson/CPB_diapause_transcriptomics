#!/bin/bash

# config_path_setup.sh - Set configuration file paths for pipeline scripts

# Log file
LOG_FILE="config_path_setup_$(date +%Y%m%d_%H%M%S).log"
echo "Starting config path setup at $(date)" > "$LOG_FILE"

# Find all shell scripts in the current directory
SCRIPTS=$(find . -maxdepth 1 -name "[0-9]*.sh" | sort)

if [ -z "$SCRIPTS" ]; then
    echo "No scripts found! Make sure you're running this from the scripts directory."
    exit 1
fi

# Get the absolute path to the current directory (where scripts are located)
SCRIPTS_DIR="$(readlink -f "$(pwd)")"
echo "Scripts directory: $SCRIPTS_DIR" | tee -a "$LOG_FILE"

# Count of modified files
MODIFIED=0

for script in $SCRIPTS; do
    # Skip this script itself
    if [[ "$(basename "$script")" == "$(basename "$0")" ]]; then
        continue
    fi
    
    # Skip the config file itself
    if [[ "$(basename "$script")" == "00.pipeline_config.sh" ]]; then
        continue
    fi
    
    echo "Processing $script..." | tee -a "$LOG_FILE"
    
    # First, check the file for problematic patterns
    CONFIG_LINE_COUNT=$(grep -c "config_file.*00.pipeline_config.sh" "$script")
    SOURCE_LINE_COUNT=$(grep -c "source.*00.pipeline_config.sh" "$script")
    
    # If multiple patterns or no patterns found, handle specially
    if [[ "$CONFIG_LINE_COUNT" -gt 1 ]] || [[ "$SOURCE_LINE_COUNT" -gt 1 ]]; then
        echo "  WARNING: Multiple config or source lines found in $script. Manual inspection required." | tee -a "$LOG_FILE"
        continue
    elif [[ "$CONFIG_LINE_COUNT" -eq 0 ]] && [[ "$SOURCE_LINE_COUNT" -eq 0 ]]; then
        echo "  No config file references found in $script. Skipping." | tee -a "$LOG_FILE"
        continue
    fi
    
    # Process the file - create a temp file with updated content
    TEMP_FILE=$(mktemp)
    
    # Read the file line by line
    while IFS= read -r line; do
        # Check for config_file pattern
        if [[ "$line" == *"config_file"*"00.pipeline_config.sh"* ]]; then
            echo "config_file=\"${SCRIPTS_DIR}/00.pipeline_config.sh\"" >> "$TEMP_FILE"
        # Check for source pattern
        elif [[ "$line" == *"source"*"00.pipeline_config.sh"* ]]; then
            echo "source \"${SCRIPTS_DIR}/00.pipeline_config.sh\"" >> "$TEMP_FILE"
        # Any other line stays as is
        else
            echo "$line" >> "$TEMP_FILE"
        fi
    done < "$script"
    
    # Replace original with temp file
    mv "$TEMP_FILE" "$script"
    chmod +x "$script"
    
    echo "  Updated $script to use path: ${SCRIPTS_DIR}/00.pipeline_config.sh" | tee -a "$LOG_FILE"
    MODIFIED=$((MODIFIED + 1))
done

echo "Config path setup completed at $(date)" | tee -a "$LOG_FILE"
echo "Modified $MODIFIED script files to use: ${SCRIPTS_DIR}/00.pipeline_config.sh" | tee -a "$LOG_FILE"
echo "See $LOG_FILE for details."