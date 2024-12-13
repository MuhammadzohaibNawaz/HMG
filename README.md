# HMG ==> GA/SA driven MDL-based DNA Compression

## Overview

This repository contains Genetic Algorithm (GA) and Simulated Annealing (SA) based approach for DNA compression, leveraging the Minimum Description Length (MDL) principle. Our method aims to efficiently compress genomic sequences by identifying optimal k-mers (patterns) that provide the best compression performance. This project addresses the challenges associated with the exponential growth of genomic data, offering a robust solution for effective data management.

## Features

- **Genetic Algorithm Optimization**: Utilizes a genetic algorithm to optimize k-mer selection, enhancing compression efficiency.
- **Simulated Annealing Optimization**: Utilizes a simulated annealing to optimize k-mer selection, enhancing compression efficiency.
- **Minimum Description Length Principle**: Applies the MDL principle to identify the most compact representation of genomic sequences.
- **Flexible Input**: Supports various genomic datasets in text and FASTA format.
- **Performance Benchmarking**: Evaluates compression ratios and time against state-of-the-art methods.

## Installation

### Prerequisites

- Java Development Kit (JDK) 8 or higher
- [Maven](https://maven.apache.org/download.cgi) (for building the project)

### Clone the Repository

### Installation ###

<pre>
git clone https://github.com/MuhammadzohaibNawaz/HMG.git
cd HMG
</pre>

### Build the Project ###

<pre>
  mvn clean install
</pre>

### Usage ###

To run the compression algorithm, modify the DNAClassification class according to your input files and parameters.
  1. **Set Dataset:** Modify the DS variable to select your dataset.
  2. **Adjust Parameters:** Change parameters such as generations, and topSubsequences (this is given as input when running the code) based on your requirements.
  3. **Run the Program:** Execute the main method to start the compression process.
<pre>
  java -cp target/HMG-1.0-SNAPSHOT.jar dna.GMG
</pre>

## Contributions 
Contributions are welcome! Please open an issue or submit a pull request for any enhancements or bug fixes.

## Acknowledgments 
1. The development of this compression algorithm was inspired by the need for efficient genomic data storage.
2. Thanks to the BioPython library for its tools and resources for handling biological data.
