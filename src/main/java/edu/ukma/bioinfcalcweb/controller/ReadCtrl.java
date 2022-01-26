package edu.ukma.bioinfcalcweb.controller;

import edu.ukma.bioinfcalcweb.service.ReadSrv;
import org.springframework.data.util.Pair;
import org.springframework.stereotype.Controller;
import org.springframework.web.bind.annotation.*;

import java.util.List;

@Controller
@RequestMapping("reads")
public class ReadCtrl {
    @PostMapping("pattern-composition")
    public String patternComposition(@RequestBody String[] patterns) {
        return ReadSrv.patternComposition(patterns);
    }

    @PostMapping("find-contigs")
    @ResponseBody
    public List<String> findContigs(@RequestBody String[] dna) {
        return ReadSrv.findContigs(dna);
    }

    // todo: check if Spring knows how to parse Pair
    @PostMapping("pattern-reconstruction")
    public String patternReconstruction(@RequestBody List<Pair<String, String>> vertices, @RequestParam int d) {
        return ReadSrv.pairedPatternReconstruction(vertices, d);
    }

}
