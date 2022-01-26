package edu.ukma.bioinfcalcweb.controller.dto;

public class StringPeptideSpectrumDto {
    private String peptide;
    private int[] spectrum;

    public String getPeptide() {
        return peptide;
    }

    public void setPeptide(String peptide) {
        this.peptide = peptide;
    }

    public int[] getSpectrum() {
        return spectrum;
    }

    public void setSpectrum(int[] spectrum) {
        this.spectrum = spectrum;
    }
}
