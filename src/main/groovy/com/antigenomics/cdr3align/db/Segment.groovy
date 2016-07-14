package com.antigenomics.cdr3align.db

class Segment {
    final String name

    Segment(String name) {
        this.name = name
    }

    boolean equals(o) {
        if (this.is(o)) return true
        if (getClass() != o.class) return false

        Segment segment = (Segment) o

        if (name != segment.name) return false

        return true
    }

    int hashCode() {
        return name.hashCode()
    }
}
