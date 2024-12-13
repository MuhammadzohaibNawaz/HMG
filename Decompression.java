package DNA;

import java.io.*;
import java.util.*;

public class Decompression {
    public static void main(String[] args) {
        String datasetName = "AgPh"; // Change this to match your dataset name
        String folderPath = "Data/geco/DS5/";
        String compressedFilePath = folderPath + "output/" + datasetName + "_compressed.bin";
        String dictionaryFilePath = folderPath + "output/" + datasetName + "_compressed.dict";
        String outputFilePath = folderPath + "output/" + datasetName + "_decompressed.txt";

        try {
            // Read dictionary file
            Map<String, Character> huffmanCodes = new HashMap<>();  // binary code -> character
            Map<Character, String> substitutions = new HashMap<>(); // special char -> original sequence
            
            System.out.println("Reading dictionary file...");
            readDictionary(dictionaryFilePath, huffmanCodes, substitutions);
            
            // Read and decode binary file
            System.out.println("Decoding binary file...");
            String decodedText = decodeBinaryFile(compressedFilePath, huffmanCodes);
            
            // Apply substitutions
            System.out.println("Applying substitutions...");
            String finalText = applySubstitutions(decodedText, substitutions);
            
            // Write output
            System.out.println("Writing decompressed file...");
            writeOutput(outputFilePath, finalText);
            
            System.out.println("Decompression completed successfully!");
            
        } catch (IOException e) {
            System.err.println("Error during decompression: " + e.getMessage());
            e.printStackTrace();
        }
    }

    private static void readDictionary(String filePath, 
                                     Map<String, Character> huffmanCodes,
                                     Map<Character, String> substitutions) throws IOException {
        try (BufferedReader reader = new BufferedReader(new FileReader(filePath))) {
            String line;
            String section = "";
            
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("===")) {
                    section = line;
                    continue;
                }
                
                if (line.trim().isEmpty()) continue;
                
                String[] parts = line.split(":");
                if (parts.length != 2) continue;
                
                if (section.contains("HUFFMAN TABLE")) {
                    // Store reversed mapping (code -> character)
                    huffmanCodes.put(parts[1], parts[0].charAt(0));
                } else if (section.contains("SUBSTITUTION TABLE")) {
                    substitutions.put(parts[0].charAt(0), parts[1]);
                }
            }
        }
    }

    private static String decodeBinaryFile(String filePath, 
                                         Map<String, Character> huffmanCodes) throws IOException {
        StringBuilder result = new StringBuilder();
        StringBuilder currentCode = new StringBuilder();
        
        try (DataInputStream in = new DataInputStream(new FileInputStream(filePath))) {
            int numSequences = in.readInt();
            
            while (in.available() > 0) {
                int currentByte = in.read();
                for (int i = 7; i >= 0; i--) {
                    boolean bit = ((currentByte >> i) & 1) == 1;
                    currentCode.append(bit ? '1' : '0');
                    
                    // Check if current code matches any Huffman code
                    if (huffmanCodes.containsKey(currentCode.toString())) {
                        result.append(huffmanCodes.get(currentCode.toString()));
                        currentCode.setLength(0);
                    }
                }
            }
        }
        return result.toString();
    }

    private static String applySubstitutions(String text, 
                                           Map<Character, String> substitutions) {
        StringBuilder result = new StringBuilder(text);
        
        // Sort substitutions by length (longest first) to handle overlapping patterns
        List<Map.Entry<Character, String>> sortedSubs = new ArrayList<>(substitutions.entrySet());
        sortedSubs.sort((a, b) -> b.getValue().length() - a.getValue().length());
        
        // Apply each substitution
        for (Map.Entry<Character, String> entry : sortedSubs) {
            int index = 0;
            while ((index = result.indexOf(String.valueOf(entry.getKey()), index)) != -1) {
                result.replace(index, index + 1, entry.getValue());
                index += entry.getValue().length();
            }
        }
        
        return result.toString();
    }

    private static void writeOutput(String filePath, String content) throws IOException {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(filePath))) {
            writer.write(content);
        }
    }
}