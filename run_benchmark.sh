#!/bin/bash
set -e

echo "Building xpress_fj..."
make

echo "Building xpress_baseline..."
(cd benchmark/baseline && make)


for timeout in 60 600
do
	echo "Running xpress_baseline $timeout second benchmarking..."
	(cd benchmark && python3 run.py --timeout $timeout baseline/xpress_baseline baseline_$timeout miplib2017benchmark)

	echo "Running xpress_fj $timeout second benchmarking..."
	(cd benchmark && python3 run.py --timeout $timeout ../xpress_fj xpress_fj_$timeout miplib2017benchmark)
done

echo "Building results report..."
(cd benchmark && python3 report.py)

echo "Benchmarking finished. See benchmark/report.html."
