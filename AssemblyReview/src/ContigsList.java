
import java.lang.reflect.Array;
import java.util.ArrayList;
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

    //List content
    //ArrayList<String> listContent = new ArrayList<String>();

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

        jList1.setModel(new javax.swing.AbstractListModel<String>() {
            String[] strings = {};
            public int getSize() { return strings.length; }
            public String getElementAt(int i) { return strings[i]; }
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

    //Set list values when fasta file is selected
    public void setList(ArrayList<StringBuffer> fastaFileContent) {
        //Instance of statisticsCalculation
        statisticsCalculation Stats = new statisticsCalculation();
        //Get arrayList of headers
        ArrayList<String> headerList = Stats.listHeaders(fastaFileContent);
        //Initialise listContent, an array of strings containing the headers
        //String[] listContent = new String[headerList.size()];
        //listContent = headerList.toArray();

        //Put values into GUI panel list
        jList1.setModel(new javax.swing.AbstractListModel<String>() {
            public int getSize() {
                return headerList.size();
            }

            public String getElementAt(int i) {
                return headerList.get(i);
            }
        });
    }


    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JList<String> jList1;
    private javax.swing.JScrollPane jScrollPane1;
    // End of variables declaration//GEN-END:variables
}
