CP.SEGREGATION.TYPES <- list("<abxcd>", "<efxeg>", "<hkxhk>", "<lmxll>", "<nnxnp>")

generateParentsMarker <- function(segType){
    genotypes <- sample(c("A", "T", "C", "G"))
    if (segType == "<abxcd>") {
        return(genotypes)
    } else if (segType == "<efxeg>") {
        parent1 <- sample(genotypes[c(1,2)])
        parent2 <- sample(genotypes[c(1,3)])
        return(c(cbind(parent1, parent2)))
    } else if (segType == "<hkxhk>") {
        parent1 <- sample(genotypes[c(1,2)])
        parent2 <- sample(genotypes[c(1,2)])
        return(c(cbind(parent1, parent2)))        
    } else if (segType == "<lmxll>") {
        parent1 <- sample(genotypes[c(1,2)])
        parent2 <- rep(genotypes[1], 2)
        return(c(cbind(parent1, parent2)))
    } else if (segType == "<nnxnp>") {
        parent1 <- rep(genotypes[1], 2)   
        parent2 <- sample(genotypes[c(1,2)])
        return(c(cbind(parent1, parent2)))
    } else if (segType == "f2 backcross") {
        parent1 <- rep(genotypes[1], 2)   
        parent2 <- c(genotypes[1], sample(genotypes[-1], 1))
        return(c(cbind(parent1, parent2)))
    } else {
        return(rep("--", 4))
    }
}

generateCPParents <- function(markersSegType, geneticPos){
    # Generate the genotypes at each marker according to the segregation type
    parentsGenotypes <- lapply(markersSegType, generateParentsMarker)
    parentsGenotypes <- matrix(unlist(parentsGenotypes), nrow = 4)
    dfGenotypes <- data.frame(marker=parentsGenotypes)
    # Construct the output dataframe
    segtypeCode <- matrix(markersSegType, nrow = 1)
    dfSegtypeCode <- data.frame(marker=segtypeCode)
    geneticPositions <- matrix(unlist(geneticPos), nrow = 1)
    dfGeneticPositions <- data.frame(marker=geneticPositions)
    cpParents <- rbind(as.matrix(dfSegtypeCode), as.matrix(dfGeneticPositions), as.matrix(dfGenotypes))
    return(cpParents)
}

generateRandomCPParents <- function(numberOfMarkers, maxGeneticDistance=130){
    markersSegType <- sample(CP.SEGREGATION.TYPES, numberOfMarkers, replace = TRUE)
    geneticPos <- sort(runif(numberOfMarkers, 0, maxGeneticDistance))
    return(generateCPParents(markersSegType, geneticPos))
}

chooseParentHaplotype <- function(haplotypes){
    # Input: a dataframe containing the haplotypes of each parent
    # Output: two dataframes, the haplotypes and the genotypes
    progenySize <- dim(haplotypes)[1]/2
    choosen1 <- sample(c(FALSE, TRUE), size = progenySize, replace = TRUE)
    choosen2 <- xor(1, choosen1)
    selected <- c(rbind(choosen1, choosen2))
    progeny <- subset(haplotypes, selected)
    s <- as.character(seq(1, dim(progeny)[1]/2))
    indivNames <- paste("indiv", s, sep = "")
    progenyHap <- data.frame(indivNames=c(rbind(indivNames, indivNames)))
    progenyHap <- cbind(progenyHap, progeny)
    progenyGeno <- aggregate(progenyHap, by=list(progenyHap$indivNames), FUN="paste0", collapse="")
    progenyGeno$indivNames <- progenyGeno$Group.1
    progenyGeno$Group.1 <- NULL
    # Sort the lines of progenyGeno
    progenyGeno <- progenyGeno[match(indivNames, progenyGeno$indivNames), ]
    # Add haplotypes information to the individual names
    hapInfo <- rep(c("_H0", "_H1"), length(indivNames))
    progenyHap$indivNames <- paste(progenyHap$indivNames, hapInfo, sep = "")
    return(list(progenyHap, progenyGeno))
}

findIndexes <- function(originalVector, positionsVector){
    # For each position "p" in position.vector
    # find the index "i" such that 
    # original.vector[i] < p and original.vector[i+1] > p
    indexesVector <- sum(originalVector < positionsVector[1])
    for (i in seq(2, length(positionsVector))){
        newIndex <- sum(originalVector < positionsVector[i])
        indexesVector <- c(indexesVector, newIndex)
    }
    return(indexesVector)
}

crossChromosomes <- function(chromosomes, nbOfCrossovers, markersPositions){
    # Choose the positions for crossovers between two given cromosomes
    # according to the markers positions
    crossedChromosomes <- chromosomes
    maxPosition <- tail(markersPositions,1)
    crosoversPos <- sort(runif(nbOfCrossovers, markersPositions[1], maxPosition))
    coIntervallLeftMarkerIndex <- findIndexes(originalVector = markersPositions, 
                                   positionsVector = crosoversPos)
    uniqueCrossoverPos <- table(coIntervallLeftMarkerIndex)
    # If unique.crossover.pos contains zero, remove it
    uniqueCrossoverPos <- uniqueCrossoverPos[as.integer(names(uniqueCrossoverPos)) > 0]
    # keep only crossovers occurring an odd number of times
    validCrossovers <- as.logical(as.integer(uniqueCrossoverPos) %% 2)
    recombLeftIndex <- as.integer(names(uniqueCrossoverPos)[validCrossovers])
    if (length(recombLeftIndex)>0) {
        endIndex <- length(markersPositions)
        for (i in seq(1, length(recombLeftIndex))){
            startIndex <- recombLeftIndex[i]+1
            temporal <- crossedChromosomes[1, startIndex:endIndex]
            crossedChromosomes[1, startIndex:endIndex] <- crossedChromosomes[2, startIndex:endIndex]
            crossedChromosomes[2, startIndex:endIndex] <- temporal
        }
    }
    return(crossedChromosomes)
}

makeCrossovers <- function(haplotypes, rho, markersPositions){
    # Imput: dataframe and a real value
    # Output: dataframe with crossed chromosomes
    crossedHaplotypes <- haplotypes
    nMax <- dim(haplotypes)[1]/2
    nbOfCrossovers <- rpois(nMax, rho)
    chromosomesWithRecomb <- subset(seq(1, nMax), nbOfCrossovers>0)
    for (i in chromosomesWithRecomb) {
        crossedHaplotypes[c(2*i-1, 2*i), ] <- crossChromosomes(haplotypes[c(2*i-1, 2*i), ], 
                                                                 nbOfCrossovers[i], markersPositions)
    }
    return(crossedHaplotypes)
}

makeGametes <- function(parentsHap, progenySize){
    gametesHap <- do.call("rbind", replicate(progenySize, parentsHap, simplify = FALSE))
    return(gametesHap)
}

generateParentsHap <- function(numOfMarkers, numOfAlleles=2){
    allelesList <- list("A", "T", "C", "G")
    p1H1 <- sample(allelesList, numOfMarkers, replace = TRUE)
    p1H2 <- sample(allelesList, numOfMarkers, replace = TRUE)
    p2H1 <- sample(allelesList, numOfMarkers, replace = TRUE)
    p2H2 <- sample(allelesList, numOfMarkers, replace = TRUE)
    parentsMatrix <- matrix(unlist(list(p1H1, p1H2, p2H1, p2H2)), ncol = numOfMarkers)
    parentsHap <- data.frame(indivNames=c("p1H1", "p1H2", "p2H1", "p2H2"),
                              marker=parentsMatrix)
    return(parentsHap)
}

simulateF1 <- function(progenySize, numOfMarkers, rho=0.5){
    parentsHap <- generateParentsHap(numOfMarkers)
    noCrossedGametes <- makeGametes(parentsHap, progenySize)
    crossedGametes <- makeCrossovers(noCrossedGametes, rho)
    progeny <- chooseParentHaplotype(crossedGametes)
    progenyHap <- progeny[1]
    progenyGeno <- progeny[2]
    return(progeny)
}

simulateCPF1 <- function(progenySize, numOfMarkers, segTypeList = NULL, 
                           possibleSegTypes = NULL, markersPositions = NULL, 
                           maxGeneticDistance=130){
    if (is.null(segTypeList)) {
        if (is.null(possibleSegTypes)) {
            possibleSegTypes <- CP.SEGREGATION.TYPES
        }
        segTypeList <- sample(possibleSegTypes, numOfMarkers, replace = TRUE)
    }
    if (is.null(markersPositions)) {
        markersPositions = sort(runif(numOfMarkers, 1, maxGeneticDistance))
    } else {
        maxGeneticDistance <- max(markersPositions)
    }
    rho <- maxGeneticDistance / 100
    # Generate the haplotypes of the parents
    parentsHap <- generateCPParents(segTypeList, markersPositions)
    # Generate pre-gametes (repeat the table of the parents as many times as individuals in the progeny)
    preGametes <- do.call("rbind", replicate(progenySize, parentsHap[-c(1,2),], simplify = FALSE))
    # Recombine pre-gametes
    recombinedGametes <- makeCrossovers(preGametes, rho, markersPositions)
    # Make the individuals by choosing from gametes
    progeny <- chooseParentHaplotype(recombinedGametes)
    progenyHap <- progeny[[1]]
    progenyGeno <- progeny[[2]]
    # Put the output in the right format
    header <- data.frame(indivNames=c("--", "--"))
    header <- cbind(header, parentsHap[c(1,2),])
    parentsHap <- as.data.frame(parentsHap[-c(1,2), ])
    parentsGeno <- parentsHap
    parentsHap$indivNames <- as.matrix(c("parent1_H0", "parent1_H1", "parent2_H0", "parent2_H1"))
    parentsHap <- parentsHap[c(numOfMarkers + 1, seq(1, numOfMarkers))]
    parentsGeno$indivNames <- as.matrix(c("parent1", "parent1", "parent2", "parent2"))
    parentsGeno <- aggregate(parentsGeno, by=list(parentsGeno$indivNames), FUN="paste0", collapse="")
    parentsGeno$Group.1 <- NULL
    parentsGeno$indivNames <- c("parent1", "parent2")
    parentsGeno <- parentsGeno[c(numOfMarkers + 1, seq(1, numOfMarkers))]
    outputHap <- rbind(header, parentsHap, progenyHap)
    outputGeno <- rbind(header, parentsGeno, progenyGeno)
    output <- list(outputHap, outputGeno)
    return(output)
}

messMrksOrder <- function(genotypes){
    # Randomly changes the order of the markers 
    # 
    # Args:
    #   genotypes: dataframe, the genotypes infromation
    #
    # Returns
    #   the same genotypes with columns in a new random order
    randomIndexes <- sample(seq(2, dim(genotypes)[2]))
    return(genotypes[, c(1, randomIndexes)])
}

shuffleLetters <- function(letters){
    temp <- sample(unlist(strsplit(letters, NULL)))
    return(paste(temp[1], temp[2], sep = ""))
}

randomOrderBases <- function(genotypes){
    # Avoid to always have the same order when converting the haplotypes into genotypes
    #
    # Args:
    #   genotypes: dataframe, the genotypes information
    #
    # Returns
    #   a dataframe with the same information as on genotypes but randomly changing the order or the letters
    temp <- as.matrix(genotypes[-c(1, 2), -c(1)])
    shuffledTemp <- mapply(shuffleLetters, temp)
    shuffledTemp <- matrix(shuffledTemp, ncol = dim(temp)[2])
    shuffledTemp <- rbind(as.matrix(genotypes[c(1, 2), -1]), shuffledTemp)
    shuffledTemp <- cbind(as.matrix(genotypes[, 1]), shuffledTemp)
    shuffledTemp <- data.frame(shuffledTemp, row.names = NULL)
    names(shuffledTemp) <- names(genotypes)
    return(shuffledTemp)
}

saveResults <- function(simResults, filenameHap = "./haplotypes.txt", 
                         filenameGeno = "./genotypes.txt"){
    hap <- simResults[[1]]
    write.table(as.matrix(hap), file = filenameHap, sep = " ", quote = FALSE, 
                col.names = FALSE, row.names = FALSE)
    gen <- simResults[[2]]
    write.table(as.matrix(gen), file = filenameGeno, sep = " ", quote = FALSE, 
                col.names = FALSE, row.names = FALSE)
}

writeRawFile <- function(rawData, fileName = "./unphased_bc_sim.raw"){
    numberOfMarkers <- dim(rawData)[1]
    numberOfIndividuals <- dim(rawData)[2] - 1
    line1 <- "data type f2 backcross"
    line2 <- paste(numberOfIndividuals, numberOfMarkers, 0, "symbols A=A H=B -=-\n")
    write.table(c(line1, line2), file = fileName, quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(rawData, file = fileName, append = TRUE, quote = FALSE, sep = "", 
                row.names = FALSE, col.names = FALSE)
}


save2iXoraFormat <- function(sim.results.genotypes, filename="./inputiXora.txt"){
    # Write the file in iXora format
    # Write the markers
    data2save <- sim.results.genotypes[-c(1,2), ]
    write.table(t(names(data2save)[-(1)]), filename, sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)
    # Write the progeny (without the last line)
    no.last.line <- as.matrix(head(data2save, -1))
    write.table(no.last.line, filename, append = TRUE, sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)
    # Write the last line
    last.line <- as.matrix(tail(data2save, 1))
    write.table(last.line, filename, append = TRUE, sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE,
                eol = "")
}

simGenotypesPhymap <- function(progSize, numOfMarkers, possibleSegTypes, maxGenDistance,
                               nbOfChromosomes, outputFileName){
    # Simulate the raw data and the physical map
    #
    # Args:
    #   progSize: integer, the size of the progeny
    #   numOfMarkers: integer, number of markers to simulate
    #   possibleSegTypes: list of character defining the segregation types in JoinMap format
    #   maxGenDistance: integer, the maximum genetic distance of on chromosome
    #   nbOfChromosomes: integer, number of chromosomes to simulate
    #   outputFileName: character, the name of the output raw file
    #
    # Returns
    #   Write two files, the raw file and the physical positions
    
    # Simulate the first chromosome
    crossingResults <- simulateCPF1(progSize, numOfMarkers, possibleSegTypes = possibleSegTypes, 
                                       maxGeneticDistance = maxGenDistance)
    
    # Convert to raw format
    genotypes <- crossingResults[[2]]
    names(genotypes)[-1] <- gsub("\\.", "-", names(genotypes)[-1])
    names(genotypes)[-1] <- paste(names(genotypes)[-1], "Chr1", sep = "") # Add the chromosome numberto the markers names
    markers <- genotypes[-c(1, 2, 3, 4), -c(1)]
    # Prepare the physical Map
    markersNames <- names(markers)
    chromosomeNumber <- rep(1, length(markersNames))
    markersPositions <- crossingResults[[2]][2, -1]
    phyMap <- rbind(markersNames, chromosomeNumber, markersPositions)
    
    # Simulate the other chromossomes
    for (i in seq(2, nbOfChromosomes)) {
        crossingResults <- simulateCPF1(progSize, numOfMarkers, possibleSegTypes = possibleSegTypes, 
                                           maxGeneticDistance = maxGenDistance)
        newGenotypes <- crossingResults[[2]]
        names(newGenotypes)[-1] <- gsub("\\.", "-", names(newGenotypes)[-1])
        names(newGenotypes)[-1] <- paste(names(newGenotypes)[-1], "Chr", i, sep = "") # Add * and the chromosome numberto the markers names
        genotypes <- cbind(genotypes, newGenotypes[, -1])
        # Prepare the physical Map
        markersNames <- names(newGenotypes)[-1]
        chromosomeNumber <- rep(i, length(markersNames))
        markersPositions <- crossingResults[[2]][2, -1]
        newPhyMap <- rbind(markersNames, chromosomeNumber, markersPositions)
        phyMap <- cbind(phyMap, newPhyMap)
    }
    
    # Convert and write to raw format
    markers <- genotypes[-c(1, 2, 3, 4), -c(1)]
    markersNames <- paste("*", names(genotypes)[-1], "\t", sep = "")
    markers.matrix <- as.matrix(markers)
    heterozygous.values <- c("AA", "TT", "CC", "GG")
    recoded <- ifelse(markers.matrix %in% heterozygous.values, "A", "B")
    recoded <- matrix(recoded, nrow = progSize)
    rawData <- data.frame(markersNames, t(recoded))
    writeRawFile(rawData, fileName=paste(outputFileName, ".raw", sep = ""))
    
    # Write the physical Map
    write.table(t(phyMap), file=paste(outputFileName, "_phyMap.txt", sep = ""), 
                row.names = FALSE, col.names = FALSE, sep = " ")
}

