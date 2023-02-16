
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.HashMap;
import javax.swing.AbstractListModel;

/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/GUIForms/JPanel.java to edit this template
 */
/**
 *
 * @author marie
 */
public class ContigsList extends javax.swing.JPanel {

    /**
     * Creates new form ContigsList
     */
    public ContigsList() {
        initComponents();
    }
    
    //Panel that will show the metrics of the selected contig
    contigsMetrics contigShowStats = new contigsMetrics();
    //Hashmap of contigs header and their corresponding sequence
    private HashMap<String, String> contigsHash = new HashMap<String, String>();

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        jScrollPane1 = new javax.swing.JScrollPane();
        jList1 = new javax.swing.JList<>();

        setBorder(javax.swing.BorderFactory.createTitledBorder(javax.swing.BorderFactory.createTitledBorder("List of contigs")));

        jList1.setBackground(new java.awt.Color(240, 240, 240));
        jList1.setModel(new javax.swing.AbstractListModel<String>() {
            String[] strings = {};
            public int getSize() { return strings.length; }
            public String getElementAt(int i) { return strings[i]; }
        });
        jList1.setSelectionMode(javax.swing.ListSelectionModel.SINGLE_SELECTION);
        jList1.addListSelectionListener(new javax.swing.event.ListSelectionListener() {
            public void valueChanged(javax.swing.event.ListSelectionEvent evt) {
                jList1ValueChanged(evt);
            }
        });
        jScrollPane1.setViewportView(jList1);

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(this);
        this.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jScrollPane1, javax.swing.GroupLayout.DEFAULT_SIZE, 273, Short.MAX_VALUE)
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jScrollPane1, javax.swing.GroupLayout.DEFAULT_SIZE, 301, Short.MAX_VALUE)
        );
    }// </editor-fold>//GEN-END:initComponents

    private void jList1ValueChanged(javax.swing.event.ListSelectionEvent evt) {//GEN-FIRST:event_jList1ValueChanged
        //Pass value to contigs metric corresponding text area, to get the metrics displayed
        contigShowStats.setSequence(contigsHash.get(jList1.getSelectedValuesList().get(0)));
        System.out.println(contigsHash.get(jList1.getSelectedValuesList()));
        System.out.println(contigsHash.get(jList1.getSelectedValuesList().get(0)));
    }//GEN-LAST:event_jList1ValueChanged

    public void setContigMetrics(contigsMetrics metrics){
        contigShowStats = metrics;
        this.contigShowStats.getSequence();
    }
    
    //Set list values when fasta file is selected
    public void setList(ArrayList<StringBuffer> fastaFileContent) {
        //Instance of statisticsCalculation
        statisticsCalculation Stats = new statisticsCalculation();
        //Get arrayList of headers, and list of sequences
        //headerList = Stats.listHeaders(fastaFileContent);
        //sequenceList = Stats.concatFasta(fastaFileContent);
        contigsHash = Stats.contigsHashMap(fastaFileContent);
        //Make and array of keys
        String[] keys = contigsHash.keySet().toArray(String[]::new);

        //Put values into GUI panel list
        jList1.setModel(new javax.swing.AbstractListModel<String>() {
            public int getSize() {
                return keys.length;
            }

            public String getElementAt(int i) {
                return keys[i];
            }
            
        });
    }


    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JList<String> jList1;
    private javax.swing.JScrollPane jScrollPane1;
    // End of variables declaration//GEN-END:variables
}
