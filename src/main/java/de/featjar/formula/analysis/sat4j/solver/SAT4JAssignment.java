/*
 * Copyright (C) 2022 Sebastian Krieter
 *
 * This file is part of formula-analysis-sat4j.
 *
 * formula-analysis-sat4j is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3.0 of the License,
 * or (at your option) any later version.
 *
 * formula-analysis-sat4j is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with formula-analysis-sat4j. If not, see <https://www.gnu.org/licenses/>.
 *
 * See <https://github.com/FeatureIDE/FeatJAR-formula-analysis-sat4j> for further information.
 */
package de.featjar.formula.analysis.sat4j.solver;

import de.featjar.formula.analysis.IAssignment;
import de.featjar.formula.analysis.bool.BooleanAssignment;

import java.util.*;

import de.featjar.formula.analysis.bool.BooleanClause;
import de.featjar.formula.analysis.bool.BooleanSolution;
import org.sat4j.core.VecInt;

/**
 * ...
 *
 * @author Sebastian Krieter
 * @author Elias Kuiter
 */
public class SAT4JAssignment implements IAssignment<Integer> {
    protected final VecInt integers;

    public SAT4JAssignment() {
        this.integers = new VecInt();
    }

    public SAT4JAssignment(BooleanAssignment assignment) {
        this.integers = new VecInt(assignment.getIntegers());
    }

    public VecInt getIntegers() {
        return integers;
    }

    public void clear() {
        integers.clear();
    }

    public void clear(int newSize) {
        integers.shrinkTo(newSize);
    }

    public void ensureSize(int size) {
        integers.ensure(size);
    }

    public int remove() {
        final int topElement = integers.get(integers.size());
        integers.pop();
        return topElement;
    }

    public void remove(int i) {
        integers.delete(i);
    }

    public void add(int var) {
        integers.push(var);
    }

    public void addAll(int[] vars) {
        integers.pushAll(new VecInt(vars));
    }

    public void addAll(BooleanAssignment assignment) {
        addAll(assignment.getIntegers());
    }

    public void replaceLast(int var) {
        integers.pop().unsafePush(var);
    }


    public void set(int index, int var) {
        integers.set(index, var);
    }

    public int[] toArray() {
        return Arrays.copyOf(integers.toArray(), integers.size());
    }

    public int[] toArray(int from) {
        return Arrays.copyOfRange(integers.toArray(), from, integers.size());
    }

    public int[] toArray(int from, int to) {
        return Arrays.copyOfRange(integers.toArray(), from, to);
    }

    public int peek() {
        return integers.get(integers.size() - 1);
    }

    public int peek(int i) {
        return integers.get(i);
    }

    public Optional<Boolean> getValue(int variable) {
        for (int i = 0; i < integers.size(); i++) {
            final int l = integers.unsafeGet(i);
            if (Math.abs(l) == variable) {
                return Optional.of(l > 0);
            }
        }
        return Optional.empty();
    }

    @Override
    public Optional<Object> getValue(Integer variable) {
        return Optional.ofNullable(getValue((int) variable).orElse(null));
    }

    public void setValue(int variable, boolean value) {
        for (int i = 0; i < integers.size(); i++) {
            final int l = integers.unsafeGet(i);
            if (Math.abs(l) == variable) {
                integers.set(i, value ? l : -l);
                return;
            }
        }
        integers.push(value ? variable : -variable);
    }

    public void removeValue(int variable) {
        for (int i = 0; i < integers.size(); i++) {
            final int l = integers.unsafeGet(i);
            if (Math.abs(l) == variable) {
                integers.delete(i);
                return;
            }
        }
    }

    @Override
    public LinkedHashMap<Integer, Object> getAll() {
        final LinkedHashMap<Integer, Object> map = new LinkedHashMap<>();
        for (int i = 0; i < integers.size(); i++) {
            final int l = integers.unsafeGet(i);
            if (l != 0) {
                map.put(Math.abs(l), l > 0);
            }
        }
        return map;
    }

    @Override
    public int size() {
        return integers.size();
    }

    @Override
    public BooleanAssignment toAssignment() {
        return new BooleanAssignment(toArray());
    }

    @Override
    public BooleanClause toClause() {
        return new BooleanClause(toArray());
    }

    @Override
    public BooleanSolution toSolution() {
        return new BooleanSolution(toArray());
    }
}
