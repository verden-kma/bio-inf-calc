package edu.ukma.bioinfcalcweb.controller;

import org.springframework.http.HttpHeaders;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.ControllerAdvice;
import org.springframework.web.bind.annotation.ExceptionHandler;
import org.springframework.web.context.request.WebRequest;
import org.springframework.web.servlet.mvc.method.annotation.ResponseEntityExceptionHandler;

@ControllerAdvice
public class BadInputExceptionHandler extends ResponseEntityExceptionHandler {

    @ExceptionHandler(Exception.class)
    public ResponseEntity<Object> handleWrongInputData(Exception unforeseen, WebRequest request) {
        return handleExceptionInternal(unforeseen, "Bad input data", new HttpHeaders(), HttpStatus.BAD_REQUEST, request);
    }
}
