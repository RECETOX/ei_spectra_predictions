#!/bin/sh

# Check if the filename argument is provided
if [ -z "$1" ]; then
  echo "Usage: $0 <filename>"
  exit 1
fi

# Replace multiple spaces with a single space
sed -E "s,[[:space:]]+, ,g" "$1" > "$1.tmp"

# Remove lines containing 'IONMODE' and 'CHARGE'
grep -v 'IONMODE\|CHARGE' "$1.tmp" > "$1.tmp1"

# Replace 'COMPOUND_NAME' with 'NAME'
sed 's/COMPOUND_NAME/NAME/g' "$1.tmp1" > "$1.tmp2"

# Move the final file to the original filename
mv "$1.tmp2" "$1"

# Remove temporary files
rm -f "$1.tmp" "$1.tmp1"