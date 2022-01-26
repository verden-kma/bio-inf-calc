package edu.ukma.bioinfcalcweb.controller;

import edu.ukma.bioinfcalcweb.controller.dto.PeptideSpectrumDto;
import edu.ukma.bioinfcalcweb.controller.dto.StringPeptideSpectrumDto;
import edu.ukma.bioinfcalcweb.service.PeptideReconstructionSrv;
import edu.ukma.bioinfcalcweb.util.Constants;
import org.springframework.stereotype.Controller;
import org.springframework.web.bind.annotation.*;

import java.util.Collection;
import java.util.List;

@Controller
@RequestMapping("peptide-reconstruction")
public class PeptideReconstructionCtrl {

    @PostMapping("theoretical-spectrum-peptide")
    @ResponseBody
    public int[] theoreticalSpectrum(@RequestBody int[] peptide,
                                     @RequestParam(required = false, defaultValue = "true") boolean isCircular) {
        return PeptideReconstructionSrv.theoreticalPeptideSpectrum(peptide, isCircular);
    }

    @PostMapping("perfect-spectrum-peptide")
    @ResponseBody
    public Collection<int[]> perfectSpectrumPeptides(@RequestBody int[] spectrum) {
        return PeptideReconstructionSrv.perfectSpectrumPeptides(spectrum);
    }

    // https://stepik.org/lesson/102/step/8?unit=8255
    @PostMapping("leaderboard-cyclopeptide")
    @ResponseBody
    public List<int[]> leaderboardCyclopeptideSequencing(@RequestBody int[] spectrum, @RequestParam int n) {
        int[] acidMasses = Constants.ACID_IDX_MASS.values().stream().mapToInt(Integer::intValue).toArray();
        return PeptideReconstructionSrv.leaderboardCyclopeptideSequencing(spectrum, acidMasses, n);
    }

    // https://stepik.org/lesson/102/step/3?unit=8255
    @PostMapping("score-peptide-string-against-spectrum")
    @ResponseBody
    public int scoreStringPeptide(@RequestBody StringPeptideSpectrumDto peptideScoring,
                                  @RequestParam(required = false, defaultValue = "true") boolean isCircular) {
        return PeptideReconstructionSrv.scorePeptide(peptideScoring.getPeptide(), peptideScoring.getSpectrum(), isCircular);
    }

    @PostMapping("score-peptide-mass-against-spectrum")
    @ResponseBody
    public int scoreMassPeptide(@RequestBody PeptideSpectrumDto peptideScoring,
                                @RequestParam(required = false, defaultValue = "true") boolean isCircular) {
        return PeptideReconstructionSrv.scorePeptide(peptideScoring.getPeptide(), peptideScoring.getSpectrum(), isCircular);
    }

    // https://stepik.org/lesson/104/step/7?unit=8257
    @PostMapping("convolution-cyclopeptide")
    @ResponseBody
    public List<int[]> convolutionCyclopeptide(@RequestBody int[] spectrum,
                                               @RequestParam int n, @RequestParam int m,
                                               @RequestParam(required = false, defaultValue = "57") int minAcidMass,
                                               @RequestParam(required = false, defaultValue = "200") int maxAcidMass) {
        return PeptideReconstructionSrv.convolutionCyclopeptideSequencing(spectrum, n, m, minAcidMass, maxAcidMass);
    }
}
