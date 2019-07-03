awk 'BEGIN {
    FS=OFS=","
} 
{
    a[$1 OFS $2 OFS $3 OFS $4 OFS $5]++ 
    b[$1 OFS $2 OFS $3 OFS $4 OFS $5]=b[$1 OFS $2 OFS $3 OFS $4 OFS $5] FILENAME "|"
    c[$1 OFS $2 OFS $3 OFS $4 OFS $5]=$0
} 
END { 
    for(i in a) 
        print b[i] "," a[i],c[i]  
}' 62403*.csv
