package edu.ukma.bioinfcalcweb.service;

import edu.ukma.bioinfcalcweb.util.Constants;
import org.springframework.data.util.Pair;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class AminoAcidSrv {
    private static final Map<Character, Character> RC_MAP = new HashMap<Character, Character>() {{
        put('A', 'U');
        put('T', 'A');
        put('C', 'G');
        put('G', 'C');
    }};

    public static String translateProtein(String rnaPattern) {
        if (rnaPattern.length() % 3 != 0) throw new IllegalArgumentException("input length must be dividable by 3");
        StringBuilder amonoAcid = new StringBuilder();
        for (int i = 0; i < rnaPattern.length(); i += 3) {
            char acid = Constants.GENETIC_CODE.get(rnaPattern.substring(i, i + 3));
            if (acid == '*') return amonoAcid.toString();
            amonoAcid.append(acid);
        }
        return amonoAcid.toString();
    }

    public static List<String> peptideEncoding(String dna, String peptide) {
        String rna = dna.replaceAll("T", "U");
        List<Pair<Integer, Integer>> dnaResPos = findPeptidePos(rna, peptide);
        StringBuilder sb = new StringBuilder(dna).reverse();
        for (int i = 0; i < sb.length(); i++) sb.setCharAt(i, RC_MAP.get(sb.charAt(i)));
        String reverseComplementRna = sb.toString();
        List<Pair<Integer, Integer>> rnaResPos = findPeptidePos(reverseComplementRna, peptide);
        List<String> res = new ArrayList<>();
        dnaResPos.forEach(posPair -> res.add(dna.substring(posPair.getFirst(), posPair.getSecond())));
        rnaResPos.forEach(posPair -> res.add(dna.substring(dna.length() - posPair.getSecond(),
                dna.length() - posPair.getFirst()).replace("U", "T")));
        return res;
    }

    private static List<Pair<Integer, Integer>> findPeptidePos(String rna, String peptide) {
        List<Pair<Integer, Integer>> res = new ArrayList<>();
        int matchIndex = 0;
        for (int i = 0; i < rna.length() - 3 * peptide.length(); i++) {
            if (Constants.GENETIC_CODE.get(rna.substring(i, i + 3)) == peptide.charAt(matchIndex)) {
                matchIndex = 1;
                int j = i + 3;
                while (true) {
                    if (matchIndex == peptide.length() - 1) {
                        res.add(Pair.of(j - (peptide.length() - 1) * 3, j + 3));
                        break;
                    }
                    if (Constants.GENETIC_CODE.get(rna.substring(j, j + 3)) == peptide.charAt(matchIndex)) {
                        matchIndex++;
                        j += 3;
                    } else break;
                }
                matchIndex = 0;
            }
        }
        return res;
    }
}
