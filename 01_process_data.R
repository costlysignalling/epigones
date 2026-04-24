# Load the RDS files
d <- readRDS("avg_docs_d2v.rds") #dimensions embeddings
dbook <- readRDS("book_metadata_1825-2000.rds") #book metadata

# Check the structure of the loaded objects
str(d)
str(dbook)

#Most represented (canonical) books
paste(dbook$author_id,dbook$shorttitle,dbook$corrected_n_of_copies)[order(dbook$corrected_n_of_copies,decreasing=T)][1:10]

#Extract surname from author_id
dbook$surname<-sub(",.*", "", dbook$author_id)

#Some IDs cannot be found in the metadata
sum(is.na(match(d$doc_id,dbook$docid)))

#This is the flipside - books without a known embedding
sum(is.na(match(dbook$docid,d$doc_id)))

#We assign those that we can
d$book<-dbook$shorttitle[match(d$doc_id,dbook$docid)]
d$name<-dbook$surname[match(d$doc_id,dbook$docid)]
d$namlong<-dbook$author_id[match(d$doc_id,dbook$docid)]

#Combine
d$ID<-paste(d$name,d$book,sep=", ")
#get rid of those, where we do not have data
d<-d[!is.na(d$name) & !is.na(d$book),]

#Only the embeddings
D<-d[,substr(names(d),1,1)=="D"]

# Compute variance for each dimension
var <- apply(D, 2, var2, na.rm = TRUE)
sum(var) #It does not sum to 1, ok...

D<-sqrt(300)*D/sqrt(sum(var)) #scale D to total variance 300, 1 per dimension

var <- apply(D, 2, var2, na.rm = TRUE)
sum(var)
