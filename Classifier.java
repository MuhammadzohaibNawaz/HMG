package DNA.classification;

import java.io.*;
import java.util.*;
import java.nio.file.Files;

public class Classifier {

    public static void main(String[] args) {
        String kmersFolderPath = "Data/classification/kmers"; 
        String datasetsFolderPath = "Data/classification/originalDatasets";
        
        File kmersFolder = new File(kmersFolderPath);
        File datasetsFolder = new File(datasetsFolderPath);
        
        // First, calculate dataset sizes
        Map<String, Long> datasetSizes = new HashMap<>();
        for (File datasetFile : datasetsFolder.listFiles()) {
            try {
                String content = new String(Files.readAllBytes(datasetFile.toPath()));
                datasetSizes.put(datasetFile.getName(), (long) content.length());
                System.out.println("Dataset " + datasetFile.getName() + 
                    " size: " + content.length() + " bases");
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        File[] kmerFiles = kmersFolder.listFiles((dir, name) -> name.endsWith("_optimal.txt"));
        if (kmerFiles != null) {
            for (File kmerFile : kmerFiles) {
                System.out.println("\nProcessing k-mer file: " + kmerFile.getName());
                Set<String> kmers = readKmersFromFile(kmerFile);
                String actualClass = kmerFile.getName().replace("_optimal.txt", "");
                System.out.println("Actual class: " + actualClass);
                
                // Store normalized frequencies for each kmer across datasets
                Map<String, Map<String, Double>> kmerDatasetFrequencies = new HashMap<>();
                
                // Process each dataset
                for (File datasetFile : datasetsFolder.listFiles()) {
                    String datasetName = datasetFile.getName();
                    long datasetSize = datasetSizes.get(datasetName);
                    System.out.println("Counting kmers in dataset: " + datasetName);
                    
                    try {
                        String content = new String(Files.readAllBytes(datasetFile.toPath()));
                        
                        for (String kmer : kmers) {
                            kmerDatasetFrequencies.putIfAbsent(kmer, new HashMap<>());
                            
                            StringBuilder modifiedContent = new StringBuilder(content);
                            int count = 0;
                            int index = 0;
                            char replacementChar = getReplacementChar(count);
                            
                            while ((index = modifiedContent.indexOf(kmer, index)) != -1) {
                                count++;
                                modifiedContent.replace(index, index + kmer.length(), 
                                    String.valueOf(replacementChar));
                                index += 1;
                            }
                            
                            // Calculate normalized frequency (percentage)
                            double frequency = (count * 100.0) / datasetSize;
                            kmerDatasetFrequencies.get(kmer).put(datasetName, frequency);
                            
                            // Debug output
                            System.out.printf("  %s: %d occurrences (%.6f%% of %s)%n", 
                                kmer, count, frequency, datasetName);
                        }
                        
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }
                
                // Classify kmers based on normalized frequencies
                int correctClassifications = 0;
                int totalKmers = kmers.size();
                
                for (String kmer : kmers) {
                    Map<String, Double> datasetFrequencies = kmerDatasetFrequencies.get(kmer);
                    String predictedClass = getPredictedClass(datasetFrequencies);
                    
                    if (predictedClass.equals(actualClass)) {
                        correctClassifications++;
                    }
                    
                    // Detailed classification output
                    System.out.println("\nKmer: " + kmer);
                    System.out.println("Frequencies across datasets:");
                    datasetFrequencies.forEach((dataset, freq) -> 
                        System.out.printf("  %s: %.6f%%%n", dataset, freq));
                    System.out.println("Predicted: " + predictedClass + 
                                     ", Actual: " + actualClass + 
                                     ", Correct: " + (predictedClass.equals(actualClass)));
                }
                
                double accuracy = (totalKmers > 0) ? 
                    (correctClassifications * 100.0 / totalKmers) : 0.0;
                System.out.printf("\nAccuracy for %s: %.2f%% (%d/%d correct)%n", 
                    kmerFile.getName(), accuracy, correctClassifications, totalKmers);
            }
        }
    }

    private static String getPredictedClass(Map<String, Double> datasetFrequencies) {
        String predictedClass = "";
        double maxFrequency = -1;
        
        for (Map.Entry<String, Double> entry : datasetFrequencies.entrySet()) {
            if (entry.getValue() > maxFrequency) {
                maxFrequency = entry.getValue();
                predictedClass = entry.getKey();
            }
        }
        
        return predictedClass;
    }

    // Read k-mers from the given file
    private static Set<String> readKmersFromFile(File kmerFile) {
        Set<String> kmers = new HashSet<>();
        try (BufferedReader br = new BufferedReader(new FileReader(kmerFile))) {
            String line;
            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (!line.isEmpty()) {
                    kmers.add(line);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return kmers;
    }

    // Count occurrences of each k-mer in a dataset file
    private static Map<String, Integer> countKmersInDataset(File datasetFile, Set<String> kmers) {
        Map<String, Integer> kmerCounts = new HashMap<>();
        try {
            String content = new String(Files.readAllBytes(datasetFile.toPath()));
            for (String kmer : kmers) {
                int count = countOccurrences(content, kmer);
                kmerCounts.put(kmer, count);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return kmerCounts;
    }

    // Count occurrences of a k-mer in a dataset string
    private static int countOccurrences(String text, String kmer) {
        int count = 0;
        int index = 0;
        while ((index = text.indexOf(kmer, index)) != -1) {
            count++;
            index += 1; // Move one position at a time to catch overlapping occurrences
        }
        return count;
    }

    // Classify the k-mers based on the frequency of their occurrence in the datasets
    private static int classifyKmers(Set<String> kmers, Map<String, Map<String, Integer>> datasetKmerCounts) {
        int correctClassifications = 0;

        for (String kmer : kmers) {
            String predictedClass = getPredictedClass(kmer, datasetKmerCounts);
            String actualClass = getActualClass(kmer);
            
            if (predictedClass.equals(actualClass)) {
                correctClassifications++;
            }
        }
        
        return correctClassifications;
    }

    // Determine the predicted class based on the highest frequency of a k-mer in a dataset
    private static String getPredictedClass(String kmer, Map<String, Map<String, Integer>> datasetKmerCounts) {
        String predictedClass = "";
        int maxCount = -1;

        for (Map.Entry<String, Map<String, Integer>> entry : datasetKmerCounts.entrySet()) {
            String datasetName = entry.getKey();
            int count = entry.getValue().getOrDefault(kmer, 0);
            
            // Debug output
            System.out.println("  " + kmer + " occurs " + count + " times in " + datasetName);
            
            if (count > maxCount) {
                maxCount = count;
                predictedClass = datasetName;
            }
        }

        return predictedClass;
    }

    // Determine the actual class for a k-mer based on the dataset
    private static String getActualClass(String kmer) {
        // Extract the actual class from the k-mer file name (assuming k-mer file names correspond to dataset names)
        return kmer.split("\\.")[0]; // this assumes your dataset files have the same name as the k-mer file
    }

    // Returns a unique character for marking counted k-mers
    private static char getReplacementChar(int count) {
        return (char) ('0' + (count % 10));
    }
}
