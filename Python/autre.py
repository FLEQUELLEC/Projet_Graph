
import gm
import geneontology as gom

# curl http://current.geneontology.org/ontology/go-basic.obo > data/go-basic.obo
print("Loading data/go-basic.obo")
coli = gom.load_OBO("Python/data/go-basic.obo")

# curl http://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/18.E_coli_MG1655.goa > data/18.E_coli_MG1655.goa
print("Loading data/18.E_coli_MG1655.goa")
gom.load_GOA(coli, "Python/data/18.E_coli_MG1655.goa")

print("Saving data/go.coli.tsv")
gm.save_delim(coli, "data/go.coli.tsv")

print("Done")
