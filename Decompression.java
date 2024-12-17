package DNA;

import java.io.*;
import java.util.*;

public class Decompression {
    public static void main(String[] args) {
        String datasetName = "sample"; // Change this to match your dataset name
        String folderPath = "Data/";
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
            int sequencesDecoded = 0;
            
            // Store all bytes in a buffer first
            byte[] allBytes = new byte[in.available()];
            in.readFully(allBytes);
            
            // Convert to a single bit string first
            StringBuilder allBits = new StringBuilder();
            for (byte b : allBytes) {
                String bits = String.format("%8s", Integer.toBinaryString(b & 0xFF)).replace(' ', '0');
                allBits.append(bits);
            }
            
            // Now process all bits sequentially without any byte boundaries
            String bitString = allBits.toString();
            int bitPos = 0;
            
            while (sequencesDecoded < numSequences && bitPos < bitString.length()) {
                currentCode.append(bitString.charAt(bitPos++));
                String code = currentCode.toString();
                
                if (huffmanCodes.containsKey(code)) {
                    char decodedChar = huffmanCodes.get(code);
                    result.append(decodedChar);
                    
                    // Debug first 20 sequences
                    if (sequencesDecoded < 20) {
                        System.out.println("Decoded #" + sequencesDecoded + ": " + decodedChar + 
                                         " from code: " + code);
                    }
                    
                    currentCode.setLength(0);
                    sequencesDecoded++;
                }
                
                // Safety check
                if (currentCode.length() > 32) {
                    System.err.println("Warning: Invalid code at bit position " + bitPos);
                    currentCode.setLength(0);
                }
            }
            
            System.out.println("Total sequences decoded: " + sequencesDecoded);
        }
        return result.toString();
    }

    private static String applySubstitutions(String text, 
                                           Map<Character, String> substitutions) {
        // Sort substitutions by length (longest first) to handle overlapping patterns
        List<Map.Entry<Character, String>> sortedSubs = new ArrayList<>(substitutions.entrySet());
        sortedSubs.sort((a, b) -> b.getValue().length() - a.getValue().length());
        
        // Process text in chunks to avoid memory issues
        int chunkSize = 1048576; // 1MB chunks
        StringBuilder finalResult = new StringBuilder();
        
        for (int start = 0; start < text.length(); start += chunkSize) {
            int end = Math.min(start + chunkSize, text.length());
            StringBuilder chunk = new StringBuilder(text.substring(start, end));
            
            // Apply substitutions to this chunk
            for (Map.Entry<Character, String> entry : sortedSubs) {
                char key = entry.getKey();
                String value = entry.getValue();
                
                int index = 0;
                while ((index = chunk.indexOf(String.valueOf(key), index)) != -1) {
                    chunk.replace(index, index + 1, value);
                    index += value.length();
                }
            }
            
            finalResult.append(chunk);
            
            // Print progress
            if (start % (chunkSize * 10) == 0) {
                System.out.printf("Processed %.2f%% of the text%n", 
                                (start * 100.0) / text.length());
            }
        }
        
        return finalResult.toString();
    }

    private static void writeOutput(String filePath, String content) throws IOException {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(filePath))) {
            writer.write(content);
        }
    }
}
