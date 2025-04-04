package DNA;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Scanner;
import java.util.Set;
import java.util.Collections;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.File;
import java.io.RandomAccessFile;

public class GACompression  {
    private static final int BITS_FOR_LENGTH = 4; // Bits to store the length of the pattern and Huffman code
    private static Map<String, String> patternToBitCode = new HashMap<>();
    private static long totalBases = 0;

    public static void main(String[] args) {
        // Define file path to read sequences
        String combination = "single_single"; // this is used for multiple variants of GA. First parameter (before_) is
                                              // for crossover and second (after_) for mutation. E.g., use
                                              // "single_single" for single point crossover and single point mutation.
        Map<String, Character> sequenceToCodeMap = new HashMap<>();
        String DS = "sample";// enter dataset name here
        List<String> chromosomeFiles = null;
        String folderPath = null;
        List<String> dnaSequences = new ArrayList<>();
        folderPath = "Data"; // path to folder
        chromosomeFiles = Arrays.asList(DS);

        boolean probBasedOccurances = true;
        for (String fileName : chromosomeFiles) {
            String filePath = folderPath + "/" + fileName;

            try {
                List<String> lines = Files.readAllLines(Paths.get(filePath));

                // Check if the file has .fasta or .fa extension
                if (fileName.endsWith(".fasta") || fileName.endsWith(".fa")) {
                    // Process FASTA format
                    StringBuilder chromosomeSequence = new StringBuilder();
                    for (int i = 0; i < lines.size(); i++) {
                        String line = lines.get(i).trim();

                        // Skip header line that starts with '>'
                        if (line.startsWith(">")) {
                            // If there is an accumulated sequence in the StringBuilder, add it to the list
                            if (chromosomeSequence.length() > 0) {
                                dnaSequences.add(chromosomeSequence.toString()); // Add sequence to list
                                chromosomeSequence.setLength(0); // Reset StringBuilder for the next sequence
                            }
                            continue; // Skip the header line
                        }
                        chromosomeSequence.append(line);
                    }

                    if (chromosomeSequence.length() > 0) {
                        dnaSequences.add(chromosomeSequence.toString());
                    }
                } else {
                    StringBuilder chromosomeSequence = new StringBuilder();
                    for (int i = 0; i < lines.size(); i++) {
                        String line = lines.get(i).trim();
                        if (i == 0 && line.startsWith(">")) {
                            continue;
                        }
                        chromosomeSequence.append(line);
                    }
                    if (chromosomeSequence.length() > 0) {
                        dnaSequences.add(chromosomeSequence.toString());
                    }
                }

            } catch (IOException e) {
                System.err.println("Error reading the file: " + fileName + " - " + e.getMessage());
                System.exit(1); // Exit the program if the file cannot be read
            }
        }
        long totalBases = 0;

        int bitsPerBase = 2; // Each base requires 2 bits

        for (String sequence : dnaSequences) {
            sequence = sequence.trim();
            long sequenceLength = sequence.length();
            totalBases += sequenceLength;
        }

        int generations = 10;
        Scanner scanner = new Scanner(System.in);
        System.out.print("Enter the number of top subsequences (m): ");
        int topSubsequences = scanner.nextInt();
        long startTime = System.currentTimeMillis(); // Get the start time
        // Set to store the unique best sequences (Ensuring no duplicates)
        String COV = null;
        String MV = null;

        Set<String> bestSequences = new HashSet<>();
        List<Integer> bestOccurrences = new ArrayList<>(); // To store corresponding occurrences

        String dna1 = null;
        String dna2 = null;

        switch (combination) {
            case "single_single":
                COV = "single";
                MV = "single";
                break;
            case "single_scramble":
                COV = "single";
                MV = "scramble";
                break;
            case "cycle_scramble":
                COV = "cycle";
                MV = "scramble";
                break;
            case "cycle_single":
                COV = "cycle";
                MV = "single";
                break;

            default:
                System.out.println("Invalid combination.");
                System.exit(0); // Exit the program
        }
        while (bestSequences.size() < topSubsequences) {

            int maxOverallOccurrences = 0;
            int secondMaxOccurrencesInGeneration = 0;
            String bestDnaString = ""; // To hold the actual DNA sequence with the most occurrences
            String secondBestDnaString = "";
            dna1 = generateRandomDNA();

            // Generate dna2 and check if it's the same as dna1
            do {
                dna2 = generateRandomDNA();
            } while (dna2.equals(dna1));

            String[] crossoverResult = null;
            String mutatedDna1 = null;
            String mutatedDna2 = null;
            // Store occurrences for each sequence along with the actual sequence
            Map<String, Integer> sequenceOccurrences = new HashMap<>();
            Map<String, String> sequenceStrings = new HashMap<>();
            sequenceOccurrences.put(dna1, countOccurrences(dnaSequences, dna1));

            sequenceStrings.put("DNA1", dna1);

            sequenceOccurrences.put(dna2, countOccurrences(dnaSequences, dna2));
            sequenceStrings.put("DNA2", dna2);
            // Loop for each generation (n generations)
            for (int gen = 1; gen <= generations; gen++) {
                if (COV.equals("single")) {
                    crossoverResult = applySinglePointCrossover(dna1, dna2);

                } else if (COV.equals("cycle")) {
                    crossoverResult = applyCycleCrossover(dna1, dna2);
                }
                if (crossoverResult[0].length() == 0 || crossoverResult[1].length() == 0) {
                    break;
                }
                if (MV.equals("single")) {
                    mutatedDna1 = applySinglePointMutation(dna1);

                    // Generate dna2 and check if it's the same as dna1
                    do {
                        mutatedDna2 = applySinglePointMutation(dna2);
                    } while (mutatedDna2.equals(mutatedDna1));
                } else if (MV.equals("scramble")) {
                    mutatedDna1 = applyScrambleMutation(dna1);
                    mutatedDna2 = applyScrambleMutation(dna2);

                }
                if (probBasedOccurances) {
                    // 50-50 chance to run either the first two lines or the last two lines
                    if (Math.random() < 0.5) {
                        // Run the first two lines
                        sequenceOccurrences.put(mutatedDna1, countOccurrences(dnaSequences, mutatedDna1));
                        sequenceStrings.put("Mutated DNA1", mutatedDna1);
                    } else {
                        // Run the last two lines
                        sequenceOccurrences.put(mutatedDna2, countOccurrences(dnaSequences, mutatedDna2));
                        sequenceStrings.put("Mutated DNA2", mutatedDna2);
                    }
                } else {
                    // If probBasedOccurrences is false, run all four lines
                    sequenceOccurrences.put(mutatedDna1, countOccurrences(dnaSequences, mutatedDna1));
                    sequenceStrings.put("Mutated DNA1", mutatedDna1);

                    sequenceOccurrences.put(mutatedDna2, countOccurrences(dnaSequences, mutatedDna2));
                    sequenceStrings.put("Mutated DNA2", mutatedDna2);
                }
                // Find the best sequence for this generation
                String bestDnaInGeneration = "";
                int maxOccurrencesInGeneration = 0;
                String secondBestDnaInGeneration = ""; // New variable for the second best
                //int secondMaxOccurrencesInGeneration = 0; // New variable for the second max occurrences

                for (Map.Entry<String, Integer> entry : sequenceOccurrences.entrySet()) {
                    if (entry.getValue() > maxOccurrencesInGeneration) {
                        // Update second best before updating the best
                        secondMaxOccurrencesInGeneration = maxOccurrencesInGeneration;
                        secondBestDnaInGeneration = bestDnaInGeneration;

                        // Update best
                        maxOccurrencesInGeneration = entry.getValue();
                        bestDnaInGeneration = entry.getKey();
                    } else if (entry.getValue() > secondMaxOccurrencesInGeneration) {
                        // Update second best only
                        secondMaxOccurrencesInGeneration = entry.getValue();
                        secondBestDnaInGeneration = entry.getKey();
                    }
                }

                // Track the overall best sequences across generations
                if (maxOccurrencesInGeneration > maxOverallOccurrences) {
                    maxOverallOccurrences = maxOccurrencesInGeneration;
                    bestDnaString = bestDnaInGeneration; // Save the best sequence
                }

                if (secondMaxOccurrencesInGeneration > 0) { // Ensure there is a valid second best
                    secondBestDnaString = secondBestDnaInGeneration; // Save the second best sequence
                }

                String firstString = null;
                String secondString = null;
                int firstCount = -1;
                int secondCount = -1;

                for (Map.Entry<String, Integer> entry : sequenceOccurrences.entrySet()) {
                    String currentString = entry.getKey();
                    int currentCount = entry.getValue();
                    if (currentCount > firstCount) {
                        // Update second place before first
                        secondString = firstString;
                        secondCount = firstCount;

                        // Update first place
                        firstString = currentString;
                        firstCount = currentCount;
                    } else if (currentCount >= secondCount) {
                        // Update second place only
                        secondString = currentString;
                        secondCount = currentCount;
                    }
                }
                dna1 = firstString;
                dna2 = secondString;
            }

            // Assign characters other than A, C, T, G for the best sequences
            // Create a StringBuilder to hold characters
            StringBuilder codeBuilder = new StringBuilder();

            // Iterate through all printable ASCII characters
            for (char c = 32; c < 127; c++) { // From ASCII 32 to 126
                // Append characters except for A, C, T, and G
                if (c != 'A' && c != 'C' && c != 'T' && c != 'G') {
                    codeBuilder.append(c);
                }
            }
            // Convert the StringBuilder to a char array
            char[] codes = codeBuilder.toString().toCharArray();

            // Ensure the best and second-best sequences are unique and do not overlap with already stored sequences
            if (maxOverallOccurrences != 0) { 
                // Handle the best sequence
                bestSequences.add(bestDnaString);
                bestOccurrences.add(maxOverallOccurrences);

                // Assign a unique code character for this bestDnaString
                sequenceToCodeMap.putIfAbsent(bestDnaString, codes[bestSequences.size() - 1]);

                // Replace this bestDnaString in dnaSequencesDummy with its code character
                for (int i = 0; i < dnaSequences.size(); i++) {
                    String sequence = dnaSequences.get(i);
                    // Replace all occurrences of bestDnaString with its code character
                    dnaSequences.set(i, sequence.replace(bestDnaString,
                            String.valueOf(sequenceToCodeMap.get(bestDnaString))));
                }
            }
            // Handle the second-best sequence if it exists
            if (secondMaxOccurrencesInGeneration > 0) { // Ensure there is a valid second best
                bestSequences.add(secondBestDnaString);
                bestOccurrences.add(secondMaxOccurrencesInGeneration);

                // Assign a unique code character for this secondBestDnaString
                sequenceToCodeMap.putIfAbsent(secondBestDnaString, codes[bestSequences.size() - 1]);

                // Replace this secondBestDnaString in dnaSequencesDummy with its code character
                for (int i = 0; i < dnaSequences.size(); i++) {
                    String sequence = dnaSequences.get(i);
                    // Replace all occurrences of secondBestDnaString with its code character
                    dnaSequences.set(i, sequence.replace(secondBestDnaString,
                            String.valueOf(sequenceToCodeMap.get(secondBestDnaString))));
                }
            }
        }

        // Modified pattern selection with dynamic length prioritization
        List<PatternInfo> sortedPatterns = new ArrayList<>();
        int i = 0;
        for (String pattern : bestSequences) {
            int frequency = bestOccurrences.get(i);
            int totalFreq = frequency;

            if (pattern.length() >= 2 && totalFreq >= 2) {
                // Calculate original size (2 bits per base)
                int originalBits = pattern.length() * 2 * totalFreq;

                // Calculate Huffman-encoded size
                int huffmanBits = estimateHuffmanSize(pattern, totalFreq);

                // Calculate dictionary overhead (pattern storage + Huffman table entry)
                int dictionaryOverhead = (pattern.length() * 2) + 8 +
                        (int) Math.ceil(Math.log(totalFreq) / Math.log(2));

                // Total compression cost
                int compressedSize = huffmanBits + dictionaryOverhead;
                int compressionBenefit = originalBits - compressedSize;

                // Combined score considering length, frequency and Huffman efficiency
                double lengthScore = pattern.length() * 2;
                double frequencyScore = Math.log(totalFreq) / Math.log(2);
                double huffmanScore = originalBits / (double) compressedSize; // compression ratio
                double combinedScore = lengthScore * frequencyScore * huffmanScore;

                // More stringent selection criteria
                if (combinedScore >= 1 && compressionBenefit > dictionaryOverhead) {
                    sortedPatterns.add(new PatternInfo(
                            pattern,
                            frequency,
                            compressionBenefit,
                            combinedScore,
                            huffmanBits));
                }
            }
            i++;
        }
        // Sort patterns considering Huffman efficiency
        Collections.sort(sortedPatterns, (a, b) -> {
            // First compare combined scores
            int scoreCompare = Double.compare(b.combinedScore, a.combinedScore);
            if (scoreCompare != 0)
                return scoreCompare;

            // If scores are equal, compare compression ratios
            double ratioA = (double) b.compressionBenefit / b.huffmanBits;
            double ratioB = (double) a.compressionBenefit / a.huffmanBits;
            return Double.compare(ratioA, ratioB);
        });

        // Process patterns in phases based on length
        List<PatternInfo> beneficialPatterns = new ArrayList<>();
        Set<String> coveredPositions = new HashSet<>();

        // First phase: longer patterns (length >= 6)
        for (PatternInfo pattern : sortedPatterns) {
            if (pattern.pattern.length() >= 6 && pattern.compressionBenefit > 100) {
                beneficialPatterns.add(pattern);
            }
        }
        for (PatternInfo pattern : sortedPatterns) {
            if (pattern.pattern.length() >= 4 && pattern.pattern.length() < 6
                    && pattern.compressionBenefit > 150) {
                beneficialPatterns.add(pattern);
            }
        }
        for (PatternInfo pattern : sortedPatterns) {
            if (pattern.pattern.length() < 4 && pattern.compressionBenefit > 200) { // Higher threshold for
                                                                                    // shorter patterns
                beneficialPatterns.add(pattern);
            }
        }

        // Limit the total number of patterns if needed
        int maxPatterns = 50; // Adjust this value based on your needs
        if (beneficialPatterns.size() > maxPatterns) {
            beneficialPatterns = beneficialPatterns.subList(0, maxPatterns);
        }

        // Generate and store bit codes for each beneficial pattern
        for (int j = 0; j < beneficialPatterns.size(); j++) {
            PatternInfo pattern = beneficialPatterns.get(j);
            patternToBitCode.put(pattern.pattern, getBitCode(j));
        }

        // Now calculate compressed size with only beneficial patterns

        long compressedBits = 0;
        int dictionaryBits = 0;

        // Calculate dictionary overhead
        for (PatternInfo pattern : beneficialPatterns) {
            String bitCode = patternToBitCode.get(pattern.pattern);
            dictionaryBits += BITS_FOR_LENGTH + (pattern.pattern.length() * bitsPerBase) + 4 + bitCode.length();
        }

        // Process sequences
        for (String sequence : dnaSequences) {
            String processedSeq = sequence;

            // Replace patterns in order of benefit
            for (PatternInfo pattern : beneficialPatterns) {
                String bitCode = patternToBitCode.get(pattern.pattern);
                int occurrences = countOccurrences(Arrays.asList(processedSeq), pattern.pattern);
                compressedBits += occurrences * bitCode.length();
                processedSeq = processedSeq.replace(pattern.pattern, "");
            }

            // Add bits for remaining bases
            compressedBits += processedSeq.length() * 2;
        }

        // Add dictionary overhead
       
        // Add file size check after compression
        File compressedFile = new File(folderPath + "/output/" + DS + "_compressed.bin");
        File dictionaryFile = new File(folderPath + "/output/" + DS + "_compressed.dict");

        // After all pattern replacements and before final statistics
        String outputPath = folderPath + "/output/" + DS + "_compressed";
        encodeAndSaveToOutput(dnaSequences, sequenceToCodeMap, outputPath);
        
        long compressedTotalBits = compressedBits + dictionaryBits;
        double bpb = (double) compressedTotalBits / totalBases;

        System.out.println("Bits per base (BPB): " + bpb);

        long endTime = System.currentTimeMillis(); // Get the start time
        double executionTimeInSeconds = (endTime - startTime) / 1000.0;
        System.out.println("Execution time: " + executionTimeInSeconds + " seconds");
    }

    // Updated PatternInfo class
    private static class PatternInfo {
        String pattern;
        int frequency;
        int compressionBenefit;
        double combinedScore;
        int huffmanBits;

        PatternInfo(String pattern, int frequency, int compressionBenefit,
                double combinedScore, int huffmanBits) {
            this.pattern = pattern;
            this.frequency = frequency;
            this.compressionBenefit = compressionBenefit;
            this.combinedScore = combinedScore;
            this.huffmanBits = huffmanBits;
        }
    }

    // Optimized bit code generation
    private static String getBitCode(int index) {
        if (index == 0)
            return "0";
        if (index == 1)
            return "1";

        int bitLength = (int) (Math.log(index) / Math.log(2)) + 1;
        StringBuilder code = new StringBuilder();

        // Create a prefix-free code
        for (int i = 0; i < bitLength - 1; i++) {
            code.append('1');
        }
        code.append('0');

        String binary = Integer.toBinaryString(index);
        code.append(binary.substring(1));

        return code.toString();
    }

    public static int countOccurrences(List<String> sequences, String target) {
        int totalOccurrences = 0;

        // Loop over all sequences
        for (String sequence : sequences) {
            sequence = sequence.trim(); // Clean up leading/trailing spaces if any
            int sequenceLength = sequence.length();
            int targetLength = target.length();

            // Use a sliding window approach to count occurrences
            for (int i = 0; i <= sequenceLength - targetLength; i++) {
                // Check if the target matches starting from the current position
                if (sequence.startsWith(target, i)) {
                    totalOccurrences++;
                }
            }
        }

        return totalOccurrences;
    }

    // Function to generate random DNA sequence of length between 2 and 6
    public static String generateRandomDNA() {
        Random rand = new Random();
        char[] nucleotides = { 'A', 'C', 'T', 'G' };

        // Generate length: size 1 with given probability, or size 2-6 otherwise
        int length = rand.nextInt(6) + 2; // Random size between 2 and 6

        // Generate the DNA sequence
        StringBuilder dna = new StringBuilder();
        for (int i = 0; i < length; i++) {
            dna.append(nucleotides[rand.nextInt(4)]); // Randomly select a nucleotide
        }
        return dna.toString();
    }

    // Function to apply single-point mutation
    public static String applySinglePointMutation(String dna) {
        if (dna.length() == 0) {
            return dna; // Return the original sequence
        }
        Random rand = new Random();
        char[] nucleotides = { 'A', 'C', 'T', 'G' };

        // Select a random position to mutate
        int mutationPosition = rand.nextInt(dna.length());
        // Find a new nucleotide different from the current one at the mutation position
        char currentBase = dna.charAt(mutationPosition);
        char newBase;
        do {
            newBase = nucleotides[rand.nextInt(4)];
        } while (newBase == currentBase); // Ensure the new base is different

        // Create the mutated DNA sequence
        StringBuilder mutatedDna = new StringBuilder(dna);
        mutatedDna.setCharAt(mutationPosition, newBase);

        return mutatedDna.toString();
    }

    // Function to apply single-point crossover between two DNA sequences
    public static String[] applySinglePointCrossover(String dna1, String dna2) {
        // System.out.println(dna1);
        //                    System.out.println(dna2);
        Random rand = new Random();

        // Ensure both sequences have the same length before crossover
        int minLength = Math.min(dna1.length(), dna2.length());
        if (minLength < 2) {
            System.out.println(
                    "Single Crossover cannot happen. One or both sequences are too short. Trying again on another sequence");
            return new String[] { "", dna2 }; // Return the original sequences
        }

        // Select a random crossover point
        int crossoverPoint = rand.nextInt(minLength);
        // Create offspring by combining the DNA from both parents at the crossover
        // point
        String offspring1 = dna1.substring(0, crossoverPoint) + dna2.substring(crossoverPoint);
        String offspring2 = dna2.substring(0, crossoverPoint) + dna1.substring(crossoverPoint);
        // Return the two offspring
        return new String[] { offspring1, offspring2 };
    }

    public static boolean isOverlapping(Set<String> bestSequences, String newSequence) {
        for (String existingSequence : bestSequences) {
            // Check if new sequence is a substring of an existing sequence, or vice versa
            if (existingSequence.contains(newSequence) || newSequence.contains(existingSequence)) {
                return true; // Overlapping detected
            }
        }
        return false; // No overlap
    }

    // Function to calculate compression size based on Shannon entropy
    public static double calculateBPBUsingShannonEntropy(Set<String> bestSequences, List<Integer> bestOccurrences) {
        // Total number of occurrences of all top subsequences
        int totalOccurrences = bestOccurrences.stream().mapToInt(Integer::intValue).sum();

        // Calculate the size of the encoded subsequences
        double totalEncodedSize = 0.0;

        // Loop through all best sequences
        for (int i = 0; i < bestSequences.size(); i++) {
            // Get the number of occurrences for the current sequence
            int occurrences = bestOccurrences.get(i);

            // Calculate the probability of the subsequence occurring
            double probability = (double) occurrences / totalOccurrences;

            // Calculate the bit length based on Shannon entropy (-log2(p))
            double bitLength = -Math.log(probability) / Math.log(2);
            totalEncodedSize += occurrences * bitLength;
        }
        return totalEncodedSize;
    }

    // Function to apply multi-point mutation
    private static char mutateNucleotide(char nucleotide) {
        Random rand = new Random();
        char[] possibleMutations = { 'A', 'T', 'C', 'G' };
        char newNucleotide;

        // Ensure the new nucleotide is different from the original one
        do {
            newNucleotide = possibleMutations[rand.nextInt(4)]; // Pick a random nucleotide
        } while (newNucleotide == nucleotide);

        return newNucleotide;
    }

    public static char getNextSymbol(Set<Character> usedSymbols) {
        // Start with printable ASCII characters that aren't DNA bases
        for (char c = 33; c < 127; c++) {
            if (c != 'A' && c != 'C' && c != 'T' && c != 'G' && !usedSymbols.contains(c)) {
                usedSymbols.add(c);
                return c;
            }
        }
        throw new RuntimeException("No more available symbols");
    }

    private static void encodeAndSaveToOutput(List<String> dnaSequences,
            Map<String, Character> sequenceToCodeMap, String outputPath) {
        try {
            // 1. Create frequency map
            Map<Character, Integer> freqMap = new HashMap<>();
            int totalSequences = 0; // Track total sequences to encode

            // First pass: count frequencies and total sequences
            for (String seq : dnaSequences) {
                for (char c : seq.toCharArray()) {
                    freqMap.merge(c, 1, Integer::sum);
                    totalSequences++;
                }
            }

            // 2. Build Huffman tree and get codes
            HuffmanTree huffman = new HuffmanTree(freqMap);
            Map<Character, String> huffmanCodes = huffman.getCodes();

            // 3. Save dictionary with clear format
            try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputPath + ".dict"))) {
                writer.write("=== SUBSTITUTION TABLE ===\n");
                for (Map.Entry<String, Character> entry : sequenceToCodeMap.entrySet()) {
                    writer.write(entry.getValue() + ":" + entry.getKey() + "\n");
                }
                writer.write("\n=== HUFFMAN TABLE ===\n");
                for (Map.Entry<Character, String> entry : huffmanCodes.entrySet()) {
                    writer.write(entry.getKey() + ":" + entry.getValue() + "\n");
                }
            }

            // 4. Save compressed data with proper bit handling
            try (DataOutputStream out = new DataOutputStream(new FileOutputStream(outputPath + ".bin"))) {
                // Write total number of sequences at the start
                out.writeInt(totalSequences);

                StringBuilder bitStream = new StringBuilder();
                int bitsWritten = 0;

                // Process each sequence
                for (String seq : dnaSequences) {
                    for (char c : seq.toCharArray()) {
                        String code = huffmanCodes.get(c);
                        if (code != null) {
                            bitStream.append(code);
                            bitsWritten += code.length();
                            StringBuilder excessBits = new StringBuilder();

                            if (bitStream.length() >= 8192) {
                                // Store the excess bits
                                int excessLength = bitStream.length() - 8192;
                                if (excessLength > 0) {
                                    // Append the excess bits to the temporary variable
                                    excessBits.append(bitStream.substring(8192));
                                }

                                // Write the first 8192 bits to the output
                                writeBitsExact(out, new StringBuilder(bitStream.substring(0, 8192)), false);

                                // Clear the bitStream for the next iteration
                                bitStream.setLength(0);
                            }

                            // At the start of the next iteration, prepend the excess bits
                            if (excessBits.length() > 0) {
                                bitStream.insert(0, excessBits);
                                excessBits.setLength(0); // Clear the excessBits after prepending
                            }
                            // Write in chunks when buffer is large enough
                            // if (bitStream.length() >= 8192) {
                            // writeBitsExact(out, bitStream, false);
                            // bitStream.setLength(0);
                            // }
                        }
                    }
                }

                // Write remaining bits with proper padding
                if (bitStream.length() > 0) {
                    writeBitsExact(out, bitStream, true);
                }
            }

            // 5. Print file sizes
            File compressedFile = new File(outputPath + ".bin");
            File dictionaryFile = new File(outputPath + ".dict");
            double compressedSizeKB = (compressedFile.length() + dictionaryFile.length()) / 1024.0;

            try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputPath + ".dict", true))) {
                writer.write("\n=== COMPRESSED SIZE ===\n");
                writer.write(String.format("Total compressed size: %.2f KB\n", compressedSizeKB));
                writer.write(String.format("Binary file: %.2f KB\n", compressedFile.length() / 1024.0));
                writer.write(String.format("Dictionary file: %.2f KB\n", dictionaryFile.length() / 1024.0));
            }

            System.out.printf("Compressed files total size: %.2f KB%n", compressedSizeKB);

        } catch (IOException e) {
            System.err.println("Error writing compressed file: " + e.getMessage());
        }
    }

    private static void writeBitsExact(DataOutputStream out, StringBuilder bitStream, boolean isLast)
            throws IOException {
        int remainingBits = bitStream.length();
        int currentByte = 0;
        int bitsInCurrentByte = 0;

        for (int i = 0; i < remainingBits; i++) {
            currentByte = (currentByte << 1) | (bitStream.charAt(i) == '1' ? 1 : 0);
            bitsInCurrentByte++;

            if (bitsInCurrentByte == 8) {
                out.write(currentByte);
                currentByte = 0;
                bitsInCurrentByte = 0;
            }
        }

        // Handle last byte with explicit padding if needed
        if (isLast && bitsInCurrentByte > 0) {
            currentByte <<= (8 - bitsInCurrentByte); // Left-align remaining bits
            out.write(currentByte);
        }
    }

    // Helper class for Huffman coding
    private static class HuffmanTree {
        private class Node implements Comparable<Node> {
            char ch;
            int freq;
            Node left, right;

            Node(char ch, int freq) {
                this.ch = ch;
                this.freq = freq;
            }

            @Override
            public int compareTo(Node other) {
                return Integer.compare(this.freq, other.freq);
            }
        }

        private final Node root;
        private final Map<Character, String> codes;

        public HuffmanTree(Map<Character, Integer> freqMap) {
            PriorityQueue<Node> pq = new PriorityQueue<>();
            freqMap.forEach((ch, freq) -> pq.offer(new Node(ch, freq)));

            while (pq.size() > 1) {
                Node left = pq.poll();
                Node right = pq.poll();
                Node parent = new Node('\0', left.freq + right.freq);
                parent.left = left;
                parent.right = right;
                pq.offer(parent);
            }
            root = pq.poll();

            codes = new HashMap<>();
            generateCodes(root, "");
        }

        private void generateCodes(Node node, String code) {
            if (node == null)
                return;
            if (node.left == null && node.right == null) {
                codes.put(node.ch, code);
                return;
            }
            generateCodes(node.left, code + "0");
            generateCodes(node.right, code + "1");
        }

        public Map<Character, String> getCodes() {
            return new HashMap<>(codes);
        }
    }

    // Add this helper method to estimate Huffman-encoded size
    private static int estimateHuffmanSize(String sequence, int frequency) {
        double probability = frequency / (double) totalBases;
        int bitsNeeded = (int) Math.ceil(-Math.log(probability) / Math.log(2));
        return bitsNeeded * frequency;
    }

    // Add this new function for cycle crossover
    private static String[] applyCycleCrossover(String dna1, String dna2) {
        int minLength = Math.min(dna1.length(), dna2.length());
        int maxLength = Math.max(dna1.length(), dna2.length());

        if (minLength < 2) {
            System.out.println("Cycle Crossover cannot happen. One or both sequences too short.");
            return new String[] { "", "" };
        }

        // Initialize offspring with the cycling portion
        char[] offspring1 = dna1.substring(0, minLength).toCharArray();
        char[] offspring2 = dna2.substring(0, minLength).toCharArray();
        boolean[] visited = new boolean[minLength];

        // Perform cycle crossover on the common length portion
        int pos = 0;
        while (!visited[pos]) {
            visited[pos] = true;
            char value = dna2.charAt(pos);
            int nextPos = dna1.substring(0, minLength).indexOf(value);

            // If value not found, break the cycle
            if (nextPos == -1)
                break;
            pos = nextPos;
        }

        // Swap unvisited positions
        for (int i = 0; i < minLength; i++) {
            if (!visited[i]) {
                char temp = offspring1[i];
                offspring1[i] = offspring2[i];
                offspring2[i] = temp;
            }
        }

        // Append remaining bases from longer sequence to both offspring
        String remainingBases = "";
        if (dna1.length() > minLength) {
            remainingBases = dna1.substring(minLength);
        } else if (dna2.length() > minLength) {
            remainingBases = dna2.substring(minLength);
        }

        return new String[] {
                new String(offspring1) + remainingBases,
                new String(offspring2) + remainingBases
        };
    }

    private static String applyScrambleMutation(String dna) {
        if (dna.length() < 2) {
            return dna; // Can't scramble a sequence of length 0 or 1
        }

        // Convert to char array for easier manipulation
        char[] dnaArray = dna.toCharArray();

        // Randomly select subset size (between 2 and length of sequence)
        Random rand = new Random();
        int subsetSize = rand.nextInt(dna.length() - 1) + 2; // At least 2 positions

        // Randomly select start position
        int startPos = rand.nextInt(dna.length() - subsetSize + 1);
        int endPos = startPos + subsetSize;

        // Scramble the selected subset
        for (int i = endPos - 1; i > startPos; i--) {
            // Pick a random index from start to i
            int j = rand.nextInt(i - startPos + 1) + startPos;

            // Swap characters at i and j
            char temp = dnaArray[i];
            dnaArray[i] = dnaArray[j];
            dnaArray[j] = temp;
        }

        return new String(dnaArray);
    }

}
