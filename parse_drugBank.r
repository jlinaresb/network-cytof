# Parsing DrugBank data
# ===

require(XML)
require(dbparser)

read_drugbank_xml_db("~/projects/networks/data/drugbank_all_full_database.xml/full database.xml")
dbparser::drugs(save_csv = T, csv_path = '~/projects/networks/data/drugbank_all_full_database.xml/')
