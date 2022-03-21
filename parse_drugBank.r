# Parsing DrugBank data
# ===

require(XML)
require(dbparser)

read_drugbank_xml_db("/mnt/netapp2/Store_uni/home/ulc/co/jlb/network-cytof/data/full database.xml")
dbparser::drugs(save_csv = T, csv_path = '/mnt/netapp2/Store_uni/home/ulc/co/jlb/network-cytof/data/')
