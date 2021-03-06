package org.opencb.opencga.core.results;

import org.opencb.commons.datastore.core.DataResult;
import org.opencb.commons.datastore.core.Event;
import org.opencb.commons.datastore.core.ObjectMap;

import java.util.ArrayList;
import java.util.List;

public class OpenCGAResult<T> extends DataResult<T> {

    private FederationNode node;

    public OpenCGAResult() {
    }

    public OpenCGAResult(int time, List<Event> events, int numResults, List<T> results, long numMatches) {
        super(time, events, numResults, results, numMatches);
    }

    public OpenCGAResult(int time, List<Event> events, int numResults, List<T> results, long numMatches, ObjectMap attributes) {
        super(time, events, numResults, results, numMatches, attributes);
    }

    public OpenCGAResult(int time, List<Event> events, long numMatches, long numInserted, long numUpdated, long numDeleted) {
        super(time, events, numMatches, numInserted, numUpdated, numDeleted);
    }

    public OpenCGAResult(int time, List<Event> events, long numMatches, long numInserted, long numUpdated, long numDeleted, 
                         ObjectMap attributes) {
        super(time, events, numMatches, numInserted, numUpdated, numDeleted, attributes);
    }

    public OpenCGAResult(int time, List<Event> events, int numResults, List<T> results, long numMatches, long numInserted, long numUpdated, 
                         long numDeleted, ObjectMap attributes) {
        super(time, events, numResults, results, numMatches, numInserted, numUpdated, numDeleted, attributes);
    }

    public OpenCGAResult(int time, List<Event> events, int numResults, List<T> results, long numMatches, long numInserted, long numUpdated,
                         long numDeleted, ObjectMap attributes, FederationNode node) {
        super(time, events, numResults, results, numMatches, numInserted, numUpdated, numDeleted, attributes);
        this.node = node;
    }

    public OpenCGAResult(DataResult<T> result) {
        this(result.getTime(), result.getEvents(), result.getNumResults(), result.getResults(), result.getNumMatches(),
                result.getNumInserted(), result.getNumUpdated(), result.getNumDeleted(), result.getAttributes());
    }

    public static OpenCGAResult empty() {
        return new OpenCGAResult<>(0, new ArrayList<>(), 0, new ArrayList<>(), 0, 0, 0, 0, new ObjectMap());
    }


    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder("OpenCGAResult{");
        sb.append("node=").append(node);
        sb.append(", time=").append(time);
        sb.append(", events=").append(events);
        sb.append(", numResults=").append(numResults);
        sb.append(", results=").append(results);
        sb.append(", resultType='").append(resultType).append('\'');
        sb.append(", numMatches=").append(numMatches);
        sb.append(", numInserted=").append(numInserted);
        sb.append(", numUpdated=").append(numUpdated);
        sb.append(", numDeleted=").append(numDeleted);
        sb.append(", attributes=").append(attributes);
        sb.append('}');
        return sb.toString();
    }

    public int getTime() {
        return time;
    }

    public OpenCGAResult<T> setTime(int time) {
        this.time = time;
        return this;
    }

    public List<Event> getEvents() {
        return events;
    }

    public OpenCGAResult<T> setEvents(List<Event> events) {
        this.events = events;
        return this;
    }

    public int getNumResults() {
        return numResults;
    }

    public OpenCGAResult<T> setNumResults(int numResults) {
        this.numResults = numResults;
        return this;
    }

    public List<T> getResults() {
        return results;
    }

    public OpenCGAResult<T> setResults(List<T> results) {
        this.results = results;
        return this;
    }

    @Deprecated
    public long getNumTotalResults() {
        return numMatches;
    }

    @Deprecated
    public OpenCGAResult<T> setNumTotalResults(long numTotalResults) {
        this.numMatches = numTotalResults;
        return this;
    }

    public long getNumMatches() {
        return numMatches;
    }

    public OpenCGAResult<T> setNumMatches(long numMatches) {
        this.numMatches = numMatches;
        return this;
    }

    public long getNumInserted() {
        return numInserted;
    }

    public OpenCGAResult<T> setNumInserted(long numInserted) {
        this.numInserted = numInserted;
        return this;
    }

    public long getNumUpdated() {
        return numUpdated;
    }

    public OpenCGAResult<T> setNumUpdated(long numUpdated) {
        this.numUpdated = numUpdated;
        return this;
    }

    public long getNumDeleted() {
        return numDeleted;
    }

    public OpenCGAResult<T> setNumDeleted(long numDeleted) {
        this.numDeleted = numDeleted;
        return this;
    }

    public String getResultType() {
        return resultType;
    }

    public OpenCGAResult<T> setResultType(String resultType) {
        this.resultType = resultType;
        return this;
    }

    public ObjectMap getAttributes() {
        return attributes;
    }

    public OpenCGAResult<T> setAttributes(ObjectMap attributes) {
        this.attributes = attributes;
        return this;
    }
    
    public FederationNode getNode() {
        return node;
    }

    public OpenCGAResult<T> setNode(FederationNode node) {
        this.node = node;
        return this;
    }
}
