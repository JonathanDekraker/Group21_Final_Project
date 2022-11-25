## Final Project



### Run Trimming

```
java -jar trimmomatic-0.39.jar SE -phred33 ../input/sars_spike_protein_raw_reads.fastq ../output.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

### Running Cleanup
```
```