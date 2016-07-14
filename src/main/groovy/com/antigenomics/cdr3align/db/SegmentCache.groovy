package com.antigenomics.cdr3align.db

import java.util.concurrent.ConcurrentHashMap
import java.util.function.Function

class SegmentCache {
    public static SegmentCache INSTANCE = new SegmentCache()

    private SegmentCache() {

    }

    private final ConcurrentHashMap<String, Segment> cache = new ConcurrentHashMap<>()

    Segment getOrCreate(String name) {
        cache.computeIfAbsent(name, new Function<String, Segment>() {
            @Override
            Segment apply(String s) {
                new Segment(name)
            }
        })
    }
}
