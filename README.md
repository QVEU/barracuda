# barracuda
 synonymous barcode designer

### Usage:
```
python Barracuda.py <sequence> [<start (1-indexed)> <maxLen>]
```

 You can run the code with:
 ```
 > qrsh
 > module load python
 ```
 then,
 ```
 > python /hpcdata/lvd_qve/QVEU_Code/barracuda/Barracuda.py ATGTTTTGA
 ```
 or you can specify position and max length (default is 30nts).
 ```
 > python /hpcdata/lvd_qve/QVEU_Code/barracuda/Barracuda.py ATGTTTTGA 3 6
 ```
