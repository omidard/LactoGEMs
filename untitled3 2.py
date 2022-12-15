



import subprocess
b = r'/Users/omidard/Desktop/reactome/groupedgprsgenes/targetseqsrxn00001.fasta'
bashCommand = "mmseqs easy-cluster b clusterRes tmp --min-seq-id 0.55 -c 0.8 --cov-mode 1"
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()