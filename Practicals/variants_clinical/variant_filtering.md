
## Filtering: Making a query database

Lets make the gemini database

```
gemini load --cores	4 -v trio.trim.vep.vcf.gz 
        -t VEP --skip-gene-tables -p recessive.ped trio.trim.vep.recessive.db
```

### Why trios?

### Variant Annotation & Impact

### Family information
