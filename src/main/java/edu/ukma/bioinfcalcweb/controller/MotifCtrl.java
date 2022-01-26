package edu.ukma.bioinfcalcweb.controller;

import edu.ukma.bioinfcalcweb.service.MotifSrv;
import org.springframework.stereotype.Controller;
import org.springframework.web.bind.annotation.*;

import java.util.List;

@Controller
@RequestMapping("motifs")
public class MotifCtrl {
    private final MotifSrv motifSrv;

    public MotifCtrl(MotifSrv motifSrv) {
        this.motifSrv = motifSrv;
    }

    @PostMapping("greedy-search")
    @ResponseBody
    public List<String> greedyMotifSearch(@RequestBody String[] dna, @RequestParam int k) {
        return MotifSrv.greedyMotifSearch(dna, k);
    }

    @PostMapping("randomized-search")
    @ResponseBody
    public List<String> randomizedMotifSearch(@RequestBody String[] dna, @RequestParam int k) {
        return motifSrv.randomizedMotifSearch(dna, k);
    }

    @PostMapping("gibbs-sampler-search")
    public List<String> gibbsSamplerSearch(@RequestBody String[] dna, @RequestParam int k, @RequestParam int n) {
        return motifSrv.gibbsSamplerSearch(dna, k, n);
    }
}
