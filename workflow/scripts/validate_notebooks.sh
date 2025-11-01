#!/bin/bash
# Script to validate notebook structure and imports
# This checks that notebooks are syntactically valid and can import the preamble

set -e

echo "Validating notebooks..."
echo "======================"

NOTEBOOK_DIR="workflow/notebooks"
WORKFLOW_DIR="workflow"

# Check if notebooks exist
if [ ! -d "$NOTEBOOK_DIR" ]; then
    echo "Error: Notebook directory not found: $NOTEBOOK_DIR"
    exit 1
fi

# Count notebooks
NOTEBOOK_COUNT=$(find "$NOTEBOOK_DIR" -name "*.ipynb" -type f | wc -l)
echo "Found $NOTEBOOK_COUNT notebook(s)"

# Validate each notebook
for notebook in "$NOTEBOOK_DIR"/*.ipynb; do
    if [ -f "$notebook" ]; then
        echo ""
        echo "Validating: $(basename "$notebook")"
        echo "---"
        
        # Check JSON structure
        if ! python3 -c "
import json
try:
    with open('$notebook') as f:
        json.load(f)
except Exception:
    exit(1)
" 2>/dev/null; then
            echo "  ✗ Invalid JSON structure"
            exit 1
        else
            echo "  ✓ Valid JSON structure"
        fi
        
        # Check for cells and structure
        VALIDATION=$(python3 -c "
import json
import sys
try:
    with open('$notebook') as f:
        data = json.load(f)
        if 'cells' not in data:
            print('Error: Missing cells key', file=sys.stderr)
            sys.exit(1)
        print(len(data['cells']))
except Exception as e:
    print(f'Error: {e}', file=sys.stderr)
    sys.exit(1)
" 2>&1)
        
        if [ $? -ne 0 ]; then
            echo "  ✗ Invalid notebook structure: $VALIDATION"
            exit 1
        fi
        echo "  ✓ Contains $VALIDATION cells"
        
        # Check for preamble import
        if grep -q "from project_utils.notebookpreamble import" "$notebook"; then
            echo "  ✓ Uses standard preamble"
        else
            echo "  ⚠ Warning: Doesn't appear to use standard preamble"
        fi
    fi
done

echo ""
echo "======================"
echo "Validation complete!"
echo ""

# Check Python files
echo "Validating Python modules..."
echo "======================"

# Check project_utils package
if [ -f "$WORKFLOW_DIR/project_utils/__init__.py" ]; then
    echo "  ✓ project_utils package exists"
else
    echo "  ✗ project_utils package not found"
    exit 1
fi

if [ -f "$WORKFLOW_DIR/project_utils/notebookpreamble.py" ]; then
    echo "  ✓ notebookpreamble module exists"
    
    # Check syntax
    if python3 -m py_compile "$WORKFLOW_DIR/project_utils/notebookpreamble.py" 2>/dev/null; then
        echo "  ✓ notebookpreamble has valid syntax"
    else
        echo "  ✗ notebookpreamble has syntax errors"
        exit 1
    fi
else
    echo "  ✗ notebookpreamble module not found"
    exit 1
fi

echo ""
echo "======================"
echo "All validations passed!"
