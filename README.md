# HCA scATAC Processing Repository

*Contact* [Caleb Lareau](mailto:clareau@broadinstitute.org)

## Step 0
```
This workflow assumes a python3 environment with several dependency packages.

Further, we assume that samtools is in the environment.
```

## Step 1
Debarcode raw reads. This step is technology specific. 

Overall, the goal to to parse out the barcode for the technology and then append it to the start of the read name. 

Some test data is provided

```
python 01_parse_scripts/jdb_sciATAC_parse_I1I2.py -a test/run1-testdata/test_R1.fastq.gz -b test/run1-testdata/test_R2.fastq.gz -i test/run1-testdata/test_I1.fastq.gz -j test/run1-testdata/test_I2.fastq.gz -o test/output/testrun1
```

## Step 2/3
Align and rebarcode with your favoriate aligner. Also, add the barcode to a sam flag instead of the read header. 

Run this for each of the split values, which you can do with a shell loop easily. 

```
sh code/02_align_reannotate.sh testrun1-split001 test/output
# calls the 03_bamReannotate.py within the shell script
```

## Step 4+
Merge all of the `.bam` files. Index the merged file. This will represent raw reads with the identified 

One option will be to run it through  `bap`, which is essential if there are bead merges to be made. 

Note on `bap`: it will append the `.bam` basename to the barcode ID, so
it's important to keep it relatively small yet informative for this 
final merged `.bam` e.g. (`20181011-Run1.bam`) will yield barcodes like (`20181011-Run1_A204P220A107P123`)