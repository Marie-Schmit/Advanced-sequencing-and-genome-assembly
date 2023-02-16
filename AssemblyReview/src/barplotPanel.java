
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.util.ArrayList;
import java.util.Collections;
import javax.swing.JPanel;

/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/GUIForms/JPanel.java to edit this template
 */
/**
 *
 * @author marie
 */
public class barplotPanel extends javax.swing.JPanel {

    /**
     * Creates new form barplotPanel
     */
    public barplotPanel() {
        initComponents();
        xCoord = new ArrayList<Integer>();
        heights = new ArrayList<Integer>();
    }

    //Bar coordinates
    private ArrayList<Integer> xCoord;
    private ArrayList<Integer> heights;
    private ArrayList<Integer> sorted_len;
    private int width;
    private int N50_index;
    private int numberBars = 0;

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        setBorder(javax.swing.BorderFactory.createTitledBorder("Barplot of contigs lenght"));

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(this);
        this.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 428, Short.MAX_VALUE)
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 258, Short.MAX_VALUE)
        );
    }// </editor-fold>//GEN-END:initComponents

    @Override
    public void paintComponent(Graphics g) {
        super.paintComponent(g);

        //N50 and other contigs have different colors
        Color N50_color = new Color(250, 128, 114);
        Color other_color = new Color(70, 130, 180);

        //Set color and rectangles
        g.setColor(other_color);

        //For every coordinate, draw a rectangle
        for (int i = 0; i < numberBars; i++) {
            if (i == N50_index) {
                g.setColor(N50_color);
                g.drawString(sorted_len.get(i).toString(), xCoord.get(i), getHeight() - heights.get(i) - 20);
            }
            else{
                g.setColor(other_color);
            }
            g.fillRect(xCoord.get(i), getHeight() - heights.get(i) - 5, width, heights.get(i));
            //Draw bar value every 5 bars
            if ((xCoord.get(i)%10 == 0)||(i==0)){
                g.drawString(sorted_len.get(i).toString(), xCoord.get(i), getHeight() - heights.get(i) - 10);
            }
        }
    }

    //Calculate each bar coordinates
    private void setCoordinates(ArrayList<Integer> sorted_len) {
        //Coordinates initialisation
        xCoord.clear();
        heights.clear();
        
        //Calculate number of bars
        if (sorted_len.size() < (this.getWidth())/1.5) {
            numberBars = sorted_len.size();
        } else {
            numberBars = (int)(this.getWidth()/1.5);
        }
        //Width calculation
        width = (int) (this.getWidth()*1.5 / (numberBars));

        //For each contig, calculate the corresponding bar x coordinate and height
        for (int i = 0; i < numberBars; i++) {
            //Set x
            double x = (i * 1.5 * width);
            xCoord.add((int)x + 10);

            //Set height
            double h = (this.getHeight()-40) * sorted_len.get(i) / sorted_len.get(0);
            heights.add((int)h);
        }

        repaint();
    }

    //Repaint the barplot when required
    public void repaintBarPlot(ArrayList<StringBuffer> fastaFileContent) {
        //Create instance of statisticsCalculation for methods
        statisticsCalculation Stats = new statisticsCalculation();
        //Calculate list of length and sort it
        sorted_len = new ArrayList<Integer>(Stats.getLength(fastaFileContent));
        //Size of list of len
        int sizeListLen = sorted_len.size();
        //Totale length of list of length
        int totalLen = sorted_len.get(sizeListLen-1); //Total length is the last element of the array list
        sorted_len.remove(sizeListLen-1);
        Collections.sort(sorted_len, Collections.reverseOrder());
        //Calculate N50
        N50_index = Stats.calculateN50(sorted_len, totalLen)[1];
        setCoordinates(sorted_len);
    }

    // Variables declaration - do not modify//GEN-BEGIN:variables
    // End of variables declaration//GEN-END:variables
}
