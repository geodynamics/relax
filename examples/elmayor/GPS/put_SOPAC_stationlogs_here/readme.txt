#eliminate dates from names of station log files
for i in *.log; do mv $i $(echo $i | sed 's/_.\{8\}//g'); done
