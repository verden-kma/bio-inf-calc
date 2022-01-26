package edu.ukma.bioinfcalcweb.controller;

import edu.ukma.bioinfcalcweb.service.AminoAcidSrv;
import org.springframework.stereotype.Controller;
import org.springframework.web.bind.annotation.PostMapping;
import org.springframework.web.bind.annotation.RequestBody;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.ResponseBody;

import java.util.List;

@Controller
@RequestMapping("amino-acid")
public class AminoAcidCtrl {

    @PostMapping("translate-protein")
    public String translateProtein(@RequestBody String rna) {
        return AminoAcidSrv.translateProtein(rna);
    }

    @PostMapping("encode-peptides")
    @ResponseBody
    public List<String> peptideEncoding(@RequestBody String[] input) {
        return AminoAcidSrv.peptideEncoding(input[0], input[1]);
    }
}
