
import java.util.ArrayList;

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
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(this);
        this.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 780, Short.MAX_VALUE)
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 300, Short.MAX_VALUE)
        );
    }// </editor-fold>//GEN-END:initComponents

    public void drawBarplot(ArrayList<StringBuffer> fastaFileContent) {

    }

    private ArrayList<Integer> getLength(ArrayList<StringBuffer> fastaFileContent) {
        //Create instance of statisticsCalculation for methods
        statisticsCalculation Stats = new statisticsCalculation();
        //Group lines by contigs, remove \n
        ArrayList<String> contigLine = Stats.concatFasta(fastaFileContent);
        //List of length
        ArrayList<Integer> list_len = new ArrayList<Integer>();

        //For each contig or header
        for (int i = 0; i < contigLine.size(); i++) {
            //Get each line as string
            String line = contigLine.get(i).toString();
            if (!line.startsWith(">")) {
                //Calculate list of length
                list_len.add(line.length());
            }
        }
        return list_len;
    }

    // Variables declaration - do not modify//GEN-BEGIN:variables
    // End of variables declaration//GEN-END:variables
}
