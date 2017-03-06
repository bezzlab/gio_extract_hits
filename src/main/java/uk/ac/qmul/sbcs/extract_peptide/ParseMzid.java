/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.qmul.sbcs.extract_peptide;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.time.LocalDateTime;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import uk.ac.ebi.jmzidml.MzIdentMLElement;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisData;
import uk.ac.ebi.jmzidml.model.mzidml.CvParam;
import uk.ac.ebi.jmzidml.model.mzidml.DBSequence;
import uk.ac.ebi.jmzidml.model.mzidml.Peptide;
import uk.ac.ebi.jmzidml.model.mzidml.PeptideEvidence;
import uk.ac.ebi.jmzidml.model.mzidml.PeptideEvidenceRef;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationItem;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationList;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationResult;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLUnmarshaller;
import java.util.HashSet;
import java.util.HashMap;
import java.util.ArrayList;
import java.util.Iterator;
import javax.xml.validation.Validator;

/**
 *
 * @author Jun Fan@qmul
 */
public class ParseMzid extends Parser {

    private static final String SCORE = "Score";
    private final String MZID_XSD = "mzIdentML1.1.0.xsd";
    private HashSet<Integer> cvTerms = new HashSet<Integer>();//the list of known CV terms, hard coded in initialize() method
    private HashSet<String> existingTerms = new HashSet<String>();
//    private HashMap<String, HashMap<String, String>> peptideCvTerms = new HashMap<String, HashMap<String, String>>();
//    private ArrayList<String> peptides = new ArrayList<String>();
//    private String fastaFile;
    private HashMap<String, String> accessions = new HashMap<String, String>();
    private HashMap<String, String> peptideSeqs = new HashMap<String, String>();
    private HashMap<String, String> pePeptide = new HashMap<String, String>();
    private HashMap<String, String> peProtein = new HashMap<String, String>();

    public ParseMzid(String input, String output) {
        this.input = input;
        this.output = output;
//        this.fastaFile = fastaFile;
        initialize();
    }

    @Override
    void parse() {
        try {
            File file = new File(input);
            System.out.println("Start at "+LocalDateTime.now());
            int mega = (int) (file.length()/1024/1024);
            System.out.println("The input file has the size of "+mega+"M");
            if(mega < 300){
                Validator validator = XMLparser.getValidator(MZID_XSD);
                boolean validFlag = XMLparser.validate(validator, input);
                if (!validFlag) {
                    System.out.println("The mzIdentML validation went wrong, program terminated");
                    System.exit(1);
                }
                System.out.println("Validation successful " + LocalDateTime.now());
            }else{
                System.out.println("It will take too long to validate a file with this size, validation skipped");
            }
//            System.exit(0);
            //main parse code is adapted from loadMzIdentML.java from X-Tracker
            MzIdentMLUnmarshaller unmarshaller = new MzIdentMLUnmarshaller(file);
            //create the mapping between protein id (e.g. DBSeq107) and accession (e.g. P12345)
            Iterator<DBSequence> dbs = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.DBSequence);
            while (dbs.hasNext()) {
                DBSequence current = dbs.next();
//line 1387 <xsd:attribute name="accession" type="xsd:string" use="required">
                accessions.put(current.getId(), current.getAccession());
            }
            System.out.println("Proteins extracted "+LocalDateTime.now());
            Iterator<Peptide> iterPeptide = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.Peptide);
            while (iterPeptide.hasNext()) {
                Peptide peptide = iterPeptide.next();
//line 1026 under PeptideType <xsd:element name="PeptideSequence" type="sequence">
                peptideSeqs.put(peptide.getId(), peptide.getPeptideSequence());
            }
            System.out.println("Peptides extracted "+LocalDateTime.now());
            Iterator<PeptideEvidence> pes = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.PeptideEvidence);
            while(pes.hasNext()){
                PeptideEvidence pe = pes.next();
//line 1276	<xsd:attribute name="dBSequence_ref" type="xsd:string" use="required">
//line 1281	<xsd:attribute name="peptide_ref" type="xsd:string" use="required">
//line 1326	<xsd:attribute name="isDecoy" type="xsd:boolean" default="false">
                //though it is optional the default is false
                if(pe.isIsDecoy()) continue;
                pePeptide.put(pe.getId(), pe.getPeptideRef());
                peProtein.put(pe.getId(), pe.getDBSequenceRef());
            }
            System.out.println("PeptideEvidence(peptide-protein relationship) extracted "+LocalDateTime.now());

            String mainScoreCV="";
            ArrayList<String> remaining = new ArrayList<String>();
            BufferedWriter out = new BufferedWriter(new FileWriter(output));
//line 367  <xsd:element name="AnalysisData" type="psi-pi:AnalysisDataType"/>
            AnalysisData analysisData = unmarshaller.unmarshal(MzIdentMLElement.AnalysisData);
//line 356  <xsd:element name="SpectrumIdentificationList" type="SpectrumIdentificationListType" maxOccurs="unbounded"/>
            List<SpectrumIdentificationList> silList = analysisData.getSpectrumIdentificationList();
            for (SpectrumIdentificationList sil : silList) {
//line 645      <xsd:element name="SpectrumIdentificationResult" type="SpectrumIdentificationResultType" maxOccurs="unbounded"/>
                List<SpectrumIdentificationResult> sirList = sil.getSpectrumIdentificationResult();
                for (SpectrumIdentificationResult sir : sirList) {
                    SpectrumIdentificationItem selected = null;
//line 823          <xsd:element name="SpectrumIdentificationItem" type="SpectrumIdentificationItemType" maxOccurs="unbounded"/>
                    List<SpectrumIdentificationItem> siiList = sir.getSpectrumIdentificationItem();
                    for (SpectrumIdentificationItem sii : siiList) {
//line 793		<xsd:attribute name="rank" type="xsd:int" use="required">
//line 798		<xsd:attribute name="passThreshold" type="xsd:boolean" use="required">
                        if (sii.isPassThreshold() && sii.getRank() == 1) {
                            selected = sii;
                            break;
                        }
                    }//end of sii list
//                  If the top SSI found 
                    if (selected != null) {
                        StringBuilder sb = new StringBuilder();
                        //If the peptide_ref is available (as it is optional)
//line 788		<xsd:attribute name="peptide_ref" type="xsd:string">
                        if (selected.getPeptideRef() != null) {
                            //PeptideEvidence captures relationship among peptide and protein
                            //multiple peptide evidences can be caused by one peptide existing in multiple proteins, e.g. isoforms
//line 760                  <xsd:element name="PeptideEvidenceRef" type="PeptideEvidenceRefType" minOccurs="1" maxOccurs="unbounded"/>
                            List<PeptideEvidenceRef> peRefList = selected.getPeptideEvidenceRef();
                            StringBuilder proSb = new StringBuilder();
                            for (PeptideEvidenceRef peRef : peRefList) {
                                if (peProtein.containsKey(peRef.getPeptideEvidenceRef())){
                                    String dbs_id = peProtein.get(peRef.getPeptideEvidenceRef()); 
                                    proSb.append(accessions.get(dbs_id));
                                    proSb.append(",");
                                }
                            }
                            if(proSb.length()>0){
                                sb.append(peptideSeqs.get(selected.getPeptideRef()));
                                sb.append("\t>");
                                proSb.deleteCharAt(proSb.length() - 1);
                                sb.append(proSb);
                                sb.append("\t\t\t");//three tabs between protein, reverse, contaminant, to score
                            }
                            //get measurements
                            HashMap<String, String> map = new HashMap<String, String>();
                            for (CvParam cv : selected.getCvParam()) {
                                if (cv.getCvRef().equals("PSI-MS") || cv.getCvRef().equals("MS")) {
                                    String[] elmts = cv.getAccession().split(":");
                                    int acc = Integer.parseInt(elmts[1]);
                                    if (cvTerms.contains(acc)) {
                                        existingTerms.add(cv.getName());
                                        map.put(cv.getName(), cv.getValue());
                                    }
                                }
                            }
                            if(mainScoreCV.length()==0){//based on the assumption that there is only one setting of used CV terms for SII
                                //find the main score for output
                                if (existingTerms.isEmpty()) {
                                    mainScoreCV = SCORE;
                                } else if (existingTerms.contains("PSM-level FDRScore")) {
                                    mainScoreCV = "PSM-level FDRScore";
                                    existingTerms.remove("PSM-level FDRScore");
                                } else if (existingTerms.contains("PSM-level combined FDRScore")) {
                                    mainScoreCV = "PSM-level combined FDRScore";
                                    existingTerms.remove("PSM-level combined FDRScore");
                                } else if (existingTerms.contains("Mascot:score")) {
                                    mainScoreCV = "Mascot:score";
                                    existingTerms.remove("Mascot:score");
                                } else {
                                    mainScoreCV = existingTerms.iterator().next();
                                    existingTerms.remove(mainScoreCV);
                                }
                                remaining.addAll(existingTerms);
                                out.append("Peptide\tProteins\tReverse\tContaminants\t");
                                out.append(mainScoreCV);
                                out.append("\tQuantitation");
                                for (String one : remaining) {
                                    out.append("\t");
                                    out.append(one);
                                }
                                out.append("\n");
                            }
                            if (mainScoreCV.equals(SCORE)){
                                sb.append("1");
                            }else{
                                if(map.containsKey(mainScoreCV)){
                                    sb.append(map.get(mainScoreCV));
                                }else{
                                    sb.append("1");
                                }
                            }
                            sb.append("\t1");//quantitation always 1 as there is no quant info in mzid
                            for (String one : remaining) {
                                sb.append("\t");
                                if (map.containsKey(one)) {
                                    sb.append(map.get(one));
                                } else {
                                    sb.append("1");
                                }
                            }
                            sb.append("\n");
                            out.append(sb);
                        }//peptide found
                    }//if sii not null
                }//end of sir list
            }
            out.flush();
            out.close();
        } catch (IOException ex) {
            Logger.getLogger(ParseMzq.class.getName()).log(Level.SEVERE, null, ex);
        }
        System.out.println("Finish at "+LocalDateTime.now());
    }

    /**
     * initialize the list of scores needed to be exported into the result tsv file
     * the list is generated by manually checking the CV library http://psidev.cvs.sourceforge.net/viewvc/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo
     * latest check: 21 Feb 2016 CV version 3.79.0 maximum CV id 1002595
     */
    private void initialize() {
//id: MS:1001143 name: search engine specific score for PSMs def: "Search engine specific peptide scores." [PSI:PI] is_a: MS:1001105 ! peptide result details
        cvTerms.add(1001143);
//id: MS:1001153 name: search engine specific score def: "Search engine specific scores." [PSI:PI] is_a: MS:1001405 ! spectrum identification result details
        cvTerms.add(1001153);
//id: MS:1001171 name: Mascot:score
        cvTerms.add(1001171);
//id: MS:1001172 name: Mascot:expectation value
        cvTerms.add(1001172);
//id: MS:1001215 name: SEQUEST:PeptideSp
        cvTerms.add(1001215);
//id: MS:1002248 name: SEQUEST:spscore
        cvTerms.add(1002248);
//id: MS:1001328 name: OMSSA:evalue
        cvTerms.add(1001328);
//id: MS:1001329 name: OMSSA:pvalue
        cvTerms.add(1001329);
//id: MS:1001419 name: SpectraST:discriminant score F  def: spectrum score
        cvTerms.add(1001419);
//id: MS:1001330 name: X\!Tandem:expect
        cvTerms.add(1001330);
//id: MS:1001390 name: Phenyx:Score
        cvTerms.add(1001390);
//id: MS:1001396 name: Phenyx:PepPvalue
        cvTerms.add(1001396);
//id: MS:1001492 name: percolator:score
        cvTerms.add(1001492);
//id: MS:1001501 name: MSFit:Mowse score
        cvTerms.add(1001501);
//id: MS:1001502 name: Sonar:Score
        cvTerms.add(1001502);
//id: MS:1001507 name: ProteinExtractor:Score
        cvTerms.add(1001507);
//id: MS:1001568 name: Scaffold:Peptide Probability
        cvTerms.add(1001568);
//id: MS:1001569 name: IdentityE Score
        cvTerms.add(1001569);
//id: MS:1001572 name: SpectrumMill:Score
        cvTerms.add(1001572);
//id: MS:1001887 name: SQID:score
        cvTerms.add(1001887);
//id: MS:1001894 name: Progenesis:confidence score
        cvTerms.add(1001894);
//id: MS:1001950 name: PEAKS:peptideScore
        cvTerms.add(1001950);
//id: MS:1001974 name: DeBunker:score
        cvTerms.add(1001974);
//id: MS:1001978 name: MSQuant:PTM-Score
        cvTerms.add(1001978);
//id: MS:1001979 name: MaxQuant:PTM Score
        cvTerms.add(1001979);
//id: MS:1001985 name: Ascore:Ascore
        cvTerms.add(1001985);
//id: MS:1002044 name: ProteinProspector:score
        cvTerms.add(1002044);
//id: MS:1002045 name: ProteinProspector:expectation value
        cvTerms.add(1002045);
//id: MS:1002052 name: MS-GF:SpecEValue
        cvTerms.add(1002052);
//id: MS:1002053 name: MS-GF:EValue
        cvTerms.add(1002053);
//id: MS:1002255 name: Comet:spscore
        cvTerms.add(1002255);
//id: MS:1002257 name: Comet:expectation value
        cvTerms.add(1002257);
//id: MS:1002262 name: Byonic:Score
        cvTerms.add(1002262);
//id: MS:1002319 name: Amanda:AmandaScore
        cvTerms.add(1002319);
//id: MS:1002338 name: Andromeda:score
        cvTerms.add(1002338);
//id: MS:1002466 name: PeptideShaker PSM score
        cvTerms.add(1002466);
//id: MS:1002545 name: xi:score
        cvTerms.add(1002545);
//id: MS:1001874 name: FDRScore
        cvTerms.add(1001874);
//id: MS:1002355 name: PSM-level FDRScore
        cvTerms.add(1002355);
//id: MS:1002356 name: PSM-level combined FDRScore
        cvTerms.add(1002356);
    }
}
