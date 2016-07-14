package com.antigenomics.cdr3align.db

import java.util.concurrent.ConcurrentHashMap

class SegmentCache {
    public static SegmentCache INSTANCE = new SegmentCache()

    private SegmentCache() {

    }

    private final ConcurrentHashMap<String, Segment> cache = new ConcurrentHashMap<>()

    Segment getOrCreate(String name) {
        cache.putIfAbsent(name, new Segment(name))
    }
}
