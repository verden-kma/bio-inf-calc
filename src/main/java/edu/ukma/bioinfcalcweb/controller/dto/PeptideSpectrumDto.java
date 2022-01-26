package edu.ukma.bioinfcalcweb.controller.dto;

public class PeptideSpectrumDto {
    private int[] peptide;
    private int[] spectrum;

    public int[] getPeptide() {
        return peptide;
    }

    public void setPeptide(int[] peptide) {
        this.peptide = peptide;
    }

    public int[] getSpectrum() {
        return spectrum;
    }

    public void setSpectrum(int[] spectrum) {
        this.spectrum = spectrum;
    }
}
