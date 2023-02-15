
import java.util.ArrayList;
import java.util.Collections;

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
    public int calculateN50(ArrayList<Integer> list_len, int totalLen) {
        ArrayList<Integer> sorted_len = list_len;
        int N50 = 0;
        int median = totalLen / 2; //Half of the genome
        int index = 0;
        int sum_len = 0;
        //Sort list by descending order
        Collections.sort(sorted_len, Collections.reverseOrder());

        //Go throught the list of length and sum up the length until half of the genome is obtained
        while (sum_len < median) {
            //Add length to sum
            sum_len += list_len.get(index);
            //Incremente index
            index++;
        }
        //N50 is the length of the contig once half of the genome is obtained
        N50 = list_len.get(index); //Index is the indice of the last contig for which the median is obtained
        return N50;
    }
    
    //GC content calculation
    //Calculate number of G and umber of C characters in one string
    private int[] numberGC(StringBuffer line) {
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
            gc = 100 * (numberG + numberC) / totalLen; //Calculate GC content
        }
        return gc;
    }
    
    
}
