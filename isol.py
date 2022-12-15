from Bio import Entrez
import pandas as pd



all_inf=pd.read_csv('/home/omidard/all_inf2.csv')

isol = []
for target in all_inf['genome_id.1']:
    acc = target
    search = Entrez.read(Entrez.esearch(db="assembly", term=acc))
    summary = Entrez.read(Entrez.esummary(db="assembly", id=search['IdList'], report="full"))
    bios = str(summary.items())
    biosa = bios.split("Biosource")[0]
    biosam = biosa.split("BioSampleId': '")[1]
    biosamp = biosam.replace("', '", "")
    biosampl = biosamp.replace('"', '')
    handle = Entrez.esummary(db='biosample', id=biosampl)
    recs = Entrez.read(handle)
    records = recs['DocumentSummarySet']['DocumentSummary'][0]
    if 'isolation source' in records['SampleData']:
        isolation_source = records['SampleData'].split('display_name="isolation source">')[1].split('</Attribute>')[0]
        isol.append(isolation_source)
    if 'isolation source' not in records['SampleData']:
        isolation_source = 'Not reported'
        isol.append(isolation_source)


all_inf['isolation_source'] = isol
all_inf.to_csv('/home/omidard/all_inf2.csv',index=False)
