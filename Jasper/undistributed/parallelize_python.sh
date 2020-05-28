#!/bin/bash
script=$1    # Which python script to execute (default: ./identify_source.py)
data=$2      # Directory where the files are stored (default: /data/astronomy-big-data/bc96e9620e41b9aba98292d37b5eec24/LoTSS_DR2_mosaic/)
file_reg=$3  # File regex (default: *.fits)
out=$4       # Output directory for files (default: /data/astronomy-big-data/bc96e9620e41b9aba98292d37b5eec24/LoTSS_DR2_writable/
batch=$5     # Number of scripts to run in parallel at a time (default: 10)
	     # The next batch is run when the prevrious one is finished, therefore the files are sorted by size, to reduce wasting time.

if [[ $1 == '' ]]; then
    script="identify_source.py"
fi

if [[ $2 == '' ]]; then
    data="/data/astronomy-big-data/bc96e9620e41b9aba98292d37b5eec24/LoTSS_DR2_mosaic"
fi

if [[ $3 == '' ]]; then
    file_reg="*.fits"
fi

if [[ $4 == '' ]]; then
    out="/data/astronomy-big-data/bc96e9620e41b9aba98292d37b5eec24/LoTSS_DR2_writable"
fi

if [[ $5 == '' ]]; then
    batch=10
fi

echo ""
echo "Script     : ${script}"
echo "Data dir   : ${data}"
echo "File regex : ${file_reg}"
echo "Out dir    : ${out}"
echo "Batches    : ${batch}"
echo ""
echo "Processing:"

count=0
for file in $(ls -dl ${data}/${file_reg} | sort -nk 5 | awk '{print $NF}'); do
    echo "${file}"
    python3 "${script}" "${file}" "${out}" &
    let count+=1
    [[ $((count%${batch})) -eq 0 ]] && wait && echo "" && echo "Processing:"
done
wait
