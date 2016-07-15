package com.antigenomics.cdr3align.db

class Segment {
    final String name

    Segment(String name) {
        this.name = name
    }

    boolean equals(o) {
        name == ((Segment) o).name
    }

    int hashCode() {
        name.hashCode()
    }
}
