
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */

/**
 *
 * @author marie
 */
public class statisticsCalculation {
    public statisticsCalculation(){
        
    }
    
    //Calculates number of sequences in the fasta file (passed as ArrayList)
    public int numberSequence(int currentNumber, String line) {
        if (line.startsWith(">")) {
            currentNumber++; //Number of sequence increments
        }
        return (currentNumber);
    }

    //Calculates lenght of one contig / scaffhold
    public int lengthSequence(int currentLength, String line) {
        if (!line.startsWith(">")) {
            currentLength += line.length(); //Length of sequence increments
        }
        return (currentLength);
    }

    //Calculate the mininum between a mim value and a string length
    public int minSequence(int currentMin, String line) {
        int lenLine = line.length();

        //If length is minimal, actualise min value
        if (lenLine < currentMin) {
            currentMin = lenLine;
        }
        return (currentMin);
    }

    //Calculate the maximum between a max value and a string length
    public int maxSequence(int currentMax, String line) {
        int lenLine = line.length();

        //If length is minimal, actualise min value
        if (lenLine > currentMax) {
            currentMax = lenLine;
        }
        return (currentMax);
    }

    //Calculates the sum of Ns between existing sum and a line that may contain Ns
    public int sumNs(int current_Ns, String line) {
        if (line.contains("N")) { //Only analyse line with letter of interest N
            for (int i = 0; i < line.length(); i++) {
                if (line.charAt(i) == ('N')) { //Count number of N
                    current_Ns++;
                }
            }
        }
        return (current_Ns);
    }

    //Calculate and return N50 from a list of contigs lengths
    public int[] calculateN50(ArrayList<Integer> list_len, int totalLen) {
        //Results contains the N50 value and its index
        int[] results = new int[2];
        
        ArrayList<Integer> sorted_len = new ArrayList<Integer>(list_len);
        int N50 = 0;
        int median = totalLen / 2; //Half of the genome
        int index = 0;
        int sum_len = 0;
        //Sort list by descending order
        Collections.sort(sorted_len, Collections.reverseOrder());

        //Go throught the list of length and sum up the length until half of the genome is obtained
        while (sum_len < median) {
            //Add length to sum
            sum_len += sorted_len.get(index);
            //Incremente index
            index++;
        }
        //N50 is the length of the contig once half of the genome is obtained
        N50 = sorted_len.get(index-1); //Index is the indice of the last contig for which the median is obtained
        //Add values to result list
        results[0] = N50;
        results[1] = index-1;
        return results;
    }
    
    //GC content calculation
    //Calculate number of G and umber of C characters in one string
    public int[] numberGC(StringBuffer line) {
        int numberG = 0;
        int numberC = 0;

        for (int i = 0; i < line.length(); i++) {
            if (!line.toString().startsWith(">")) { //Only consider sequences
                //Count number of g or G in the line
                if (line.charAt(i) == ('G')) { //If character is a G
                    numberG++;
                } //Count number of C and c in the line
                else if (line.charAt(i) == ('C')) {
                    numberC++;
                }
            }
        }
        int counts[] = {numberG, numberC};
        return counts; //Return a list of number of G and number of C in the considered line
    }

    //Calculate G/C content of all the sequences of a fasta file
    public double getGC(ArrayList<StringBuffer> fileContent, int totalLen) {
        int numberG = 0; //Number of G in all the sequences
        int numberC = 0; //Number of C in all the sequences
        double gc; //Value of GC content

        for (int i = 0; i < fileContent.size(); i++) {
            //Get a list of the number of G and C in the line if the line is a sequence
            int counts[] = new int[2];
            counts = numberGC(fileContent.get(i));

            //Actualise number of G and C in all the sequences
            numberG += counts[0];
            numberC += counts[1];
        }

        //Division by 0 if number of C is null
        if (totalLen == 0) {
            gc = 0;
            //Throw exception
            throw new IllegalStateException("Division by 0. The length of the fasta file sequence is null. File might be empty, please try with a new file.");
        } else {
            gc = 100 * (numberG + numberC) / (double)totalLen; //Calculate GC content
        }
        return gc;
    }
    
    //Save each contig (lines between two header) as a single string in an array list.
    public ArrayList<String> concatFasta(ArrayList<StringBuffer> fileContent) {
        ArrayList<String> sequenceContent = new ArrayList<String>();
        String newLine = "";

        //For each line of the file
        for (int i = 0; i < fileContent.size(); i++) {
            //Get each line of the file and convert to string
            String line = fileContent.get(i).toString();
            //Remove line return
            line = line.replace("\n", "");

            if (!line.startsWith(">")) { //Concatenate lines that are not headers
                newLine += line;
            } else { //When header met, save concatenation and header in ArrayList
                if(newLine != "")
                    sequenceContent.add(new String(newLine));
                sequenceContent.add(new String(line));
                newLine = "";
            }
        }
        //Add last sequence
        if(newLine != "")
                    sequenceContent.add(new String(newLine));
        
        return sequenceContent;
    }
    
    //Store every header in an array list
    public ArrayList<String> listHeaders(ArrayList<StringBuffer> fileContent) {
        ArrayList<String> headers = new ArrayList<String>();
        //For each line of the file, store the one starting with ">"
        for (int i = 0; i < fileContent.size(); i++) {
            //Get each line of the file and convert to string
            String line = fileContent.get(i).toString();
            if(line.startsWith(">")){
                headers.add(line);
            }
        }
        return headers;
    }
    
    //Get a list of length of the entered fasta file
    public ArrayList<Integer> getLength(ArrayList<StringBuffer> fastaFileContent) {
        //Create instance of statisticsCalculation for methods
        statisticsCalculation Stats = new statisticsCalculation();
        //Group lines by contigs, remove \n
        ArrayList<String> contigLine = Stats.concatFasta(fastaFileContent);
        //List of length
        ArrayList<Integer> list_len = new ArrayList<Integer>();
        //Total length
        int totalLen = 0;

        //For each contig or header
        for (int i = 0; i < contigLine.size(); i++) {
            //Get each line as string
            String line = contigLine.get(i).toString();
            if (!line.startsWith(">")) {
                //Calculate list of length
                list_len.add(line.length());
                //Caclulate total length
                totalLen = lengthSequence(totalLen, line);
            }
        }
        list_len.add(totalLen); //Last value of the list is the total length
        return list_len;
    }
    
    
    //Get a hashmap of contigs. Each key is a contig header
    public HashMap<String, String> contigsHashMap(ArrayList<StringBuffer> fastaFileContent) {
        //Hashmap initialisation
        HashMap<String, String> contigsHash = new HashMap<String, String>();
        
        //For each contigs, add hader as key and sequence as value
        String newLine = ""; //Sequence value
        String key = "";

        //For each line of the file
        for (int i = 0; i < fastaFileContent.size(); i++) {
            //Get each line of the file and convert to string
            String line = fastaFileContent.get(i).toString();

            if (line.startsWith(">")) { //Header
                //New header met: end of previous contig, add contig and key value to HashMap
                contigsHash.put(new String(key), new String(newLine));
                newLine = "";
                //Create new key for current header
                key = line;
            } else { //Add lines between two headers to the same sequence
                newLine += line;
            }
        }
        return contigsHash;
    }
    
    //Get a string of contigs and a string of contigs header
    public String[][] contigsArrays(ArrayList<StringBuffer> fastaFileContent) {
        int lenFile = fastaFileContent.size();
        //Strings initialisation
        String[][] contigsAll = new String[2][lenFile];
        String
        
        //For each contigs, add header as key and sequence as value
        String newLine = ""; //Sequence value
        String key = "";

        //For each line of the file
        for (int i = 0; i < fastaFileContent.size(); i++) {
            //Get each line of the file and convert to string
            String line = fastaFileContent.get(i).toString();

            if (line.startsWith(">")) { //Header
                //New header met: end of previous contig, add contig and key value to HashMap
                contigsHash.put(new String(key), new String(newLine));
                newLine = "";
                //Create new key for current header
                key = line;
            } else { //Add lines between two headers to the same sequence
                newLine += line;
            }
        }
        return contigsHash;
    }
    
}
