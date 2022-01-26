package edu.ukma.bioinfcalcweb.util;

import com.google.common.collect.BiMap;
import com.google.common.collect.HashBiMap;

import java.io.IOException;
import java.net.URISyntaxException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

public class Constants {
    public final static Map<String, Character> GENETIC_CODE = new HashMap<>();

    public static final Map<Character, Integer> AMINO_ACID_MASS = new HashMap<>();

    public static final int MAX_MASS;

    public static final BiMap<Integer, Integer> ACID_IDX_MASS = HashBiMap.create();

    static {
        // populate GENETIC_CODE
        try {
            Files.lines(Paths.get(Constants.class.getClassLoader().getResource("static/rna_codon_table.txt").toURI()))
                    .forEach(line -> {
                        String[] codonAcid = line.split("\\s");
                        assert (codonAcid.length == 2);
                        GENETIC_CODE.put(codonAcid[0], codonAcid[1].charAt(0));
                    });
        } catch (IOException | URISyntaxException e) {
            e.printStackTrace();
            System.exit(-2);
        }

        // populate AMINO_ACID_MASS
        try {
            Files.lines(Paths.get(Constants.class.getClassLoader().getResource("static/int_mass_table.txt").toURI()))
                    .forEach(line -> {
                        String[] aam = line.split("\\s");
                        assert (aam.length == 2);
                        AMINO_ACID_MASS.put(aam[0].charAt(0), Integer.parseInt(aam[1]));
                    });
        } catch (IOException | URISyntaxException e) {
            e.printStackTrace();
            System.exit(-2);
        }
        // init ACID_IDX_MASS
        for (Integer acidMass : new HashSet<>(AMINO_ACID_MASS.values()))
            ACID_IDX_MASS.put(ACID_IDX_MASS.size(), acidMass);

        MAX_MASS = AMINO_ACID_MASS.values().stream().max(Integer::compareTo).orElse(-1);
    }
}