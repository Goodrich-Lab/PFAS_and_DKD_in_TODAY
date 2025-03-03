

// Diabetic sinvle nodes alone:
MATCH (d:Disease {commonName: "Diabetic Nephropathy"}) RETURN d
MATCH (c:Chemical {xrefDTXSID: "DTXSID8031863"}) RETURN c
MATCH (g:Gene {geneSymbol: "COLGALT1"}) RETURN g
MATCH (d:Disease {commonName: "Diabetic Nephropathy"}) RETURN d
MATCH (da:Disease {commonName: "Kidney Failure, Chronic"})


// This returns a graph with all CHEMICALINCREASESEXPRESSION associations with PFNA
MATCH (c:Chemical)-[:CHEMICALINCREASESEXPRESSION|CHEMICALDECREASESEXPRESSION]-(Gene)
WHERE c.xrefDTXSID = 'DTXSID8031863' 
RETURN c, Gene



// This is a great query!! find all shortest paths beetween genes linked to PFNA with the outcome, limited to 2 steps
MATCH (c:Chemical)-[:CHEMICALINCREASESEXPRESSION|CHEMICALDECREASESEXPRESSION]-(gene:Gene)
WHERE c.xrefDTXSID = 'DTXSID8031863' 
WITH gene, c
MATCH (d:Disease {commonName: "Diabetic Nephropathy"}), p=allShortestPaths((d)-[*..2]-(gene))
RETURN c, p;

// Look at all 3 way hops between the genes of interest and the outcome
UNWIND ["PSENEN", "SGK1", "APOA1", "GPNMB", "TGIF2LX", "ATG3", "GDF11", "MSTN", "IL271", "FGF8", "TMEM87B", "RNF114", "COLGALT1"] AS mygenelist  
MATCH (gene:Gene {geneSymbol: mygenelist})
WITH gene
MATCH (d:Disease {commonName: "Diabetic Nephropathy"}), p=allShortestPaths((d)-[*..3]-(gene))  
RETURN p LIMIT 50;


// Look at all 3 way hops between the genes of interest and PFNA
UNWIND ["PSENEN", "SGK1", "APOA1", "GPNMB", "TGIF2LX", "ATG3", "GDF11", "MSTN", "IL271", "FGF8", "TMEM87B", "RNF114", "COLGALT1"] AS mygenelist  
MATCH (gene:Gene {geneSymbol: mygenelist})
WITH gene
MATCH (c:Chemical {xrefDTXSID: "DTXSID8031863"}), p=allShortestPaths((gene)-[*..3]-(c))  
RETURN p LIMIT 50;


// Look at all 3 way hops between the exposure of outcomes of interest
UNWIND ["Diabetic Nephropathy", "Kidney Failure, Chronic"] AS myoutcomes  
MATCH (d:Disease {commonName: myoutcomes})
WITH d
MATCH (c:Chemical {xrefDTXSID: "DTXSID8031863"}), p=allShortestPaths((d)-[*..3]-(c))  
RETURN p LIMIT 50;


// Limit to only gene nodes
UNWIND ["PSENEN", "SGK1", "APOA1", "GPNMB", "TGIF2LX", "ATG3", "GDF11", "MSTN", "IL271", "FGF8", "TMEM87B", "RNF114", "COLGALT1"] AS mygenelist  
MATCH (gene:Gene {geneSymbol: mygenelist})
WITH gene
MATCH (d:Disease {commonName: "Diabetic Nephropathy"}), p=allShortestPaths((d)-[*..3]-(gene))  
WHERE all(node in nodes(p)[1..-1] WHERE node:Gene )
RETURN p LIMIT 50;




// This is a great query!! find all shortest paths beetween genes linked to PFNA with the 2 outcomes, limited to 2 steps
MATCH (c:Chemical)-[:CHEMICALINCREASESEXPRESSION|CHEMICALDECREASESEXPRESSION]-(gene:Gene)
WHERE c.xrefDTXSID = 'DTXSID8031863' 
WITH c, collect(gene) AS genes
MATCH (dn:Disease {commonName: "Diabetic Nephropathy"})
MATCH (da:Disease {commonName: "Kidney Failure, Chronic"})
WITH c, genes, dn, da
UNWIND genes AS gene
OPTIONAL MATCH p1=allShortestPaths((dn)-[*..2]-(gene))
OPTIONAL MATCH p2=allShortestPaths((da)-[*..2]-(gene))
RETURN c, gene, p1, p2;




// This is an equivelant graph, except it doesn't include nodes associatied with only PFNA and nothing else
MATCH (c:Chemical)-[:CHEMICALINCREASESEXPRESSION|CHEMICALDECREASESEXPRESSION]-(gene:Gene)
WHERE c.xrefDTXSID = 'DTXSID8031863' 
WITH gene, c
MATCH (d:Disease {commonName: "Diabetic Nephropathy"}), p=allShortestPaths((d)-[*..2]-(gene))
RETURN c, p;



















// Playing around with subgraphall 

MATCH (c:Chemical {xrefDTXSID: "DTXSID8031863"})
CALL apoc.path.subgraphAll(c, {
	relationshipFilter: "CHEMICALINCREASESEXPRESSION|CHEMICALDECREASESEXPRESSION",
    minLevel: 1,
    maxLevel: 2
})
YIELD nodes, relationships
RETURN c, nodes, relationships;






//Final Query, return as path
MATCH (d:Disease {commonName: "Diabetic Nephropathy"})
MATCH (end:Gene)
WHERE end.geneSymbol IN ["PSENEN", "SGK1", "APOA1", "GPNMB", "TGIF2LX", "ATG3", "GDF11", "MSTN", "IL271", "FGF8", "TMEM87B", "RNF114", "COLGALT1"]
CALL apoc.path.spanningTree(d, {
	relationshipFilter: "<GENEASSOCIATESWITHDISEASE|GENEINTERACTSWITHGENE",
    minLevel: 1,
    maxLevel: 5, 
    endNodes: end
})
YIELD path
RETURN path
UNION
MATCH (c:Chemical {xrefDTXSID: "DTXSID8031863"})
MATCH (d:Disease {commonName: "Diabetic Nephropathy"})
MATCH (end:Gene)
WHERE end.geneSymbol IN ["PSENEN", "SGK1", "APOA1", "GPNMB", "TGIF2LX", "ATG3", "GDF11", "MSTN", "IL271", "FGF8", "TMEM87B", "RNF114", "COLGALT1"]
CALL apoc.path.spanningTree(c, {
	relationshipFilter: "CHEMICALINCREASESEXPRESSION>|CHEMICALDECREASESEXPRESSION>|GENEINTERACTSWITHGENE",
    minLevel: 1,
    maxLevel: 7, 
    endNodes: end
})
YIELD path
RETURN path

// Final query return as nodes
MATCH (d:Disease {commonName: "Diabetic Nephropathy"})
MATCH (end:Gene)
WHERE end.geneSymbol IN ["PSENEN", "SGK1", "APOA1", "GPNMB", "TGIF2LX", "ATG3", "GDF11", "MSTN", "IL271", "FGF8", "TMEM87B", "RNF114", "COLGALT1"]
CALL apoc.path.spanningTree(d, {
	relationshipFilter: "<GENEASSOCIATESWITHDISEASE|GENEINTERACTSWITHGENE",
    minLevel: 1,
    maxLevel: 5, 
    endNodes: end
})
YIELD path
WITH nodes(path) as n, relationships(path) as p
RETURN n, p
UNION
MATCH (c:Chemical {xrefDTXSID: "DTXSID8031863"})
MATCH (d:Disease {commonName: "Diabetic Nephropathy"})
MATCH (end:Gene)
WHERE end.geneSymbol IN ["PSENEN", "SGK1", "APOA1", "GPNMB", "TGIF2LX", "ATG3", "GDF11", "MSTN", "IL271", "FGF8", "TMEM87B", "RNF114", "COLGALT1"]
CALL apoc.path.spanningTree(c, {
	relationshipFilter: "CHEMICALINCREASESEXPRESSION>|CHEMICALDECREASESEXPRESSION>|GENEINTERACTSWITHGENE",
    minLevel: 1,
    maxLevel: 7, 
    endNodes: end
})
YIELD path
WITH nodes(path) as n, relationships(path) as p
RETURN n, p


// Final query with gephi streaming
MATCH (d:Disease {commonName: "Diabetic Nephropathy"})
MATCH (end:Gene)
WHERE end.geneSymbol IN ["PSENEN", "SGK1", "APOA1", "GPNMB", "TGIF2LX", "ATG3", "GDF11", "MSTN", "IL271", "FGF8", "TMEM87B", "RNF114", "COLGALT1"]
CALL apoc.path.spanningTree(d, {
	relationshipFilter: "<GENEASSOCIATESWITHDISEASE|GENEINTERACTSWITHGENE",
    minLevel: 1,
    maxLevel: 5, 
    endNodes: end
})
YIELD path
RETURN d, path
UNION
MATCH (c:Chemical {xrefDTXSID: "DTXSID8031863"})
MATCH (d:Disease {commonName: "Diabetic Nephropathy"})
MATCH (end:Gene)
WHERE end.geneSymbol IN ["PSENEN", "SGK1", "APOA1", "GPNMB", "TGIF2LX", "ATG3", "GDF11", "MSTN", "IL271", "FGF8", "TMEM87B", "RNF114", "COLGALT1"]
CALL apoc.path.spanningTree(c, {
	relationshipFilter: "CHEMICALINCREASESEXPRESSION>|CHEMICALDECREASESEXPRESSION>|GENEINTERACTSWITHGENE",
    minLevel: 1,
    maxLevel: 7, 
    endNodes: end
})
YIELD path
call apoc.gephi.add(null,'test_workspace', path) yield nodes, relationships, time
RETURN d, path



// Query with ge
MATCH (d:Disease {commonName: "Diabetic Nephropathy"})
MATCH (end:Gene)
WHERE end.geneSymbol IN ["PSENEN", "SGK1", "APOA1", "GPNMB", "TGIF2LX", "ATG3", "GDF11", "MSTN", "IL271", "FGF8", "TMEM87B", "RNF114", "COLGALT1"]
CALL apoc.path.spanningTree(d, {
	relationshipFilter: "<GENEASSOCIATESWITHDISEASE|GENEINTERACTSWITHGENE",
    minLevel: 1,
    maxLevel: 5, 
    endNodes: end
})
YIELD path
RETURN path
UNION
MATCH (c:Chemical {xrefDTXSID: "DTXSID8031863"})
MATCH (d:Disease {commonName: "Diabetic Nephropathy"})
MATCH (end:Gene)
WHERE end.geneSymbol IN ["PSENEN", "SGK1", "APOA1", "GPNMB", "TGIF2LX", "ATG3", "GDF11", "MSTN", "IL271", "FGF8", "TMEM87B", "RNF114", "COLGALT1"]
CALL apoc.path.spanningTree(c, {
	relationshipFilter: "CHEMICALINCREASESEXPRESSION>|CHEMICALDECREASESEXPRESSION>|GENEINTERACTSWITHGENE",
    minLevel: 1,
    maxLevel: 7, 
    endNodes: end
})
YIELD path
RETURN path

CALL apoc.export.graphml.data(path, [], "movies-l.graphml", {})
YIELD file, source, format, nodes, relationships, properties, time, rows, batchSize, batches, done, data
RETURN file, source, format, nodes, relationships, properties, time, rows, batchSize, batches, done, data










call apoc.export.graphml.data([], path, 'queryRelationship.graphml',
  {useTypes: true, source: {id: 'name'}, label: {id: 'age'}})
