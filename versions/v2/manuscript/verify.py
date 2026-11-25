import json
from fractions import Fraction
from itertools import combinations_with_replacement, product
from math import prod


def determinant(matrix):
    work = [[Fraction(value) for value in row] for row in matrix]
    result = Fraction(1)
    for column in range(len(work)):
        pivot = next((row for row in range(column, len(work)) if work[row][column]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            result = -result
        value = work[column][column]
        result *= value
        for row in range(column + 1, len(work)):
            factor = work[row][column] / value
            for other in range(column + 1, len(work)):
                work[row][other] -= factor * work[column][other]
    return result


def rank(matrix):
    work = [[Fraction(value) for value in row] for row in matrix]
    pivot_row = 0
    for column in range(len(work[0]) if work else 0):
        pivot = next((row for row in range(pivot_row, len(work)) if work[row][column]), None)
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        for row in range(pivot_row + 1, len(work)):
            factor = work[row][column] / work[pivot_row][column]
            for other in range(column + 1, len(work[0])):
                work[row][other] -= factor * work[pivot_row][other]
            work[row][column] = Fraction(0)
        pivot_row += 1
        if pivot_row == len(work):
            break
    return pivot_row


def inertia(matrix):
    work = [[Fraction(value) for value in row] for row in matrix]
    positive = zero = negative = 0
    while work:
        size = len(work)
        pivot = next((index for index in range(size) if work[index][index]), None)
        if pivot is not None:
            order = [pivot] + [index for index in range(size) if index != pivot]
            work = [[work[row][column] for column in order] for row in order]
            value = work[0][0]
            positive += value > 0
            negative += value < 0
            work = [[work[row][column] - work[row][0] * work[0][column] / value
                     for column in range(1, size)] for row in range(1, size)]
            continue
        pair = next(((row, column) for row in range(size)
                     for column in range(row + 1, size) if work[row][column]), None)
        if pair is None:
            zero += size
            break
        first, second = pair
        order = [first, second] + [index for index in range(size) if index not in pair]
        work = [[work[row][column] for column in order] for row in order]
        value = work[0][1]
        positive += 1
        negative += 1
        work = [[work[row][column] -
                 (work[row][0] * work[1][column] + work[row][1] * work[0][column]) / value
                 for column in range(2, size)] for row in range(2, size)]
    return [positive, zero, negative]


def source_blocks(sizes):
    factors = [3 * size - 4 for size in sizes]

    def phi(indices):
        return prod(factors[index] for index in indices)

    def psi(indices):
        return sum(sizes[index] * phi([other for other in indices if other != index])
                   for index in indices)

    indices = list(range(len(sizes)))
    denominator = psi(indices)
    assert denominator
    diagonal = []
    within = []
    between = []
    for index, size in enumerate(sizes):
        remaining = [other for other in indices if other != index]
        omitted_psi = psi(remaining)
        omitted_theta = phi(remaining) + omitted_psi
        diagonal.append(Fraction((size - 1) * omitted_theta +
                                 2 * (size - 3) * omitted_psi, 4 * denominator))
        within.append(Fraction(-omitted_theta - 2 * omitted_psi, 4 * denominator))
        between.append([Fraction(phi([other for other in remaining if other != column]),
                                 denominator) if column != index else Fraction(0)
                        for column in indices])
    return diagonal, within, between


def full_matrices(sizes):
    diagonal, within, between = source_blocks(sizes)
    vertices = [(part, vertex) for part, size in enumerate(sizes) for vertex in range(size)]
    distance = []
    laplacian = []
    for first in vertices:
        distance_row = []
        laplacian_row = []
        for second in vertices:
            if first == second:
                distance_row.append(0)
                laplacian_row.append(diagonal[first[0]])
            elif first[0] == second[0]:
                distance_row.append(4)
                laplacian_row.append(within[first[0]])
            else:
                distance_row.append(1)
                laplacian_row.append(between[first[0]][second[0]])
        distance.append(distance_row)
        laplacian.append(laplacian_row)
    return distance, laplacian


def shift(matrix, scalar, weights=None):
    weights = weights or [1] * len(matrix)
    return [[value - (scalar * weights[row] if row == column else 0)
             for column, value in enumerate(values)] for row, values in enumerate(matrix)]


def charpoly(matrix):
    size = len(matrix)
    work = [[Fraction(row == column) for column in range(size)] for row in range(size)]
    coefficients = [Fraction(1)]
    for degree in range(1, size + 1):
        product = [[sum(matrix[row][inner] * work[inner][column] for inner in range(size))
                    for column in range(size)] for row in range(size)]
        coefficient = -sum(product[index][index] for index in range(size)) / degree
        coefficients.append(coefficient)
        work = [[value + (coefficient if row == column else 0)
                 for column, value in enumerate(values)] for row, values in enumerate(product)]
    return coefficients


def polynomial_from_roots(roots):
    coefficients = [Fraction(1)]
    for root in roots:
        following = [Fraction(0)] * (len(coefficients) + 1)
        for index, coefficient in enumerate(coefficients):
            following[index] += coefficient
            following[index + 1] -= coefficient * root
        coefficients = following
    return coefficients


def matrix_product(first, second):
    return [[sum(first[row][inner] * second[inner][column] for inner in range(len(second)))
             for column in range(len(second[0]))] for row in range(len(first))]


def polynomial_product(first, second):
    result = [Fraction(0)] * (len(first) + len(second) - 1)
    for first_index, first_value in enumerate(first):
        for second_index, second_value in enumerate(second):
            result[first_index + second_index] += first_value * second_value
    return result


def secular_coefficients(diagonal, weights, scale):
    result = polynomial_from_roots(diagonal)
    for index, weight in enumerate(weights):
        omitted = polynomial_from_roots([value for other, value in enumerate(diagonal) if other != index])
        for degree, coefficient in enumerate(omitted, 1):
            result[degree] -= scale * weight * coefficient
    return result


def polynomial_remainder(dividend, divisor):
    work = [Fraction(value) for value in dividend]
    while work and len(work) >= len(divisor):
        factor = work[0] / divisor[0]
        for index, value in enumerate(divisor):
            work[index] -= factor * value
        while work and work[0] == 0:
            work.pop(0)
    return work


def sturm_sequence(coefficients):
    degree = len(coefficients) - 1
    derivative = [value * (degree - index) for index, value in enumerate(coefficients[:-1])]
    sequence = [coefficients, derivative]
    while sequence[-1]:
        remainder = polynomial_remainder(sequence[-2], sequence[-1])
        if not remainder:
            break
        sequence.append([-value for value in remainder])
    return sequence


def sign_variations(sequence, endpoint):
    signs = []
    for coefficients in sequence:
        if endpoint == '-infinity':
            value = coefficients[0] * (-1) ** (len(coefficients) - 1)
        elif endpoint == '+infinity':
            value = coefficients[0]
        else:
            value = Fraction(0)
            for coefficient in coefficients:
                value = value * endpoint + coefficient
        if value:
            signs.append(1 if value > 0 else -1)
    return sum(first != second for first, second in zip(signs, signs[1:]))


def graph_extension_checks():
    distance_cases = inverse_cases = laplacian_polynomials = zero_cofactor_cases = 0
    negative_spectra = []
    for part_count in range(2, 9):
        for sizes in combinations_with_replacement(range(1, 9), part_count):
            total = sum(sizes)
            if total > 24:
                continue
            nonsingletons = sum(size >= 2 for size in sizes)
            singletons = part_count - nonsingletons
            factors = [3 * size - 4 for size in sizes]
            ratio = sum((Fraction(size, factor) for size, factor in zip(sizes, factors)), Fraction(0))
            distance_form = [[sizes[row] * sizes[column] +
                              (sizes[row] * factors[row] if row == column else 0)
                              for column in range(part_count)] for row in range(part_count)]
            distance_quotient = [[Fraction(value, sizes[row]) for value in values]
                                 for row, values in enumerate(distance_form)]
            quotient_inertia = inertia(distance_form)
            expected = ([nonsingletons + 1, 0, total - nonsingletons - 1] if ratio < -1 else
                        [nonsingletons, 1, total - nonsingletons - 1] if ratio == -1 else
                        [nonsingletons, 0, total - nonsingletons])
            assert [quotient_inertia[0], quotient_inertia[1],
                    quotient_inertia[2] + total - part_count] == expected
            assert charpoly(distance_quotient) == secular_coefficients(factors, sizes, Fraction(1))
            assert rank(shift(distance_form, Fraction(-1), sizes)) == part_count - max(singletons - 1, 0)
            distance_cases += 1
            if ratio != -1:
                inverse_quotient = [[(Fraction(1, factors[row]) if row == column else 0) -
                                     Fraction(sizes[column], factors[row] * factors[column]) / (1 + ratio)
                                     for column in range(part_count)] for row in range(part_count)]
                assert matrix_product(distance_quotient, inverse_quotient) == [
                    [int(row == column) for column in range(part_count)] for row in range(part_count)]
                inverse_cases += 1
            if ratio == 0:
                assert Fraction(part_count, 4) < singletons <= Fraction(part_count, 2)
                if 2 * singletons == part_count:
                    assert all(size <= 2 for size in sizes)
                zero_cofactor_cases += 1
                continue
            diagonal, within, between = source_blocks(sizes)
            laplacian_quotient = [[diagonal[row] + (sizes[row] - 1) * within[row] if row == column
                                   else sizes[column] * between[row][column]
                                   for column in range(part_count)] for row in range(part_count)]
            polynomial = secular_coefficients([Fraction(-1, factor) for factor in factors],
                                               [Fraction(size, factor * factor) for size, factor in zip(sizes, factors)],
                                               1 / ratio)
            assert charpoly(laplacian_quotient) == polynomial
            identity_product = matrix_product(laplacian_quotient, distance_quotient)
            assert all(identity_product[row][column] + int(row == column) ==
                       Fraction(sizes[column], factors[row]) / ratio
                       for row in range(part_count) for column in range(part_count))
            laplacian_polynomials += 1
            if ratio != -1:
                continue
            assert Fraction(part_count + 3, 4) < singletons <= Fraction(part_count + 1, 2)
            if 2 * singletons == part_count + 1:
                assert all(size <= 2 for size in sizes)
            distance, laplacian = full_matrices(sizes)
            assert rank(distance) == total - 1
            assert inertia(distance) == expected
            base_determinant = (-4) ** (total - part_count) * prod(factors)
            assert determinant([values[1:] for values in distance[1:]]) == Fraction(-base_determinant, factors[0] ** 2)
            kernel_vector = [Fraction(1, factors[part]) for part, size in enumerate(sizes) for vertex in range(size)]
            assert all(sum(value * coordinate for value, coordinate in zip(values, kernel_vector)) == 0
                       for values in distance)
            assert charpoly(laplacian) == polynomial_product(
                polynomial_from_roots([Fraction(1, 4)] * (total - part_count)), polynomial)
            grouped_sizes = sorted(set(sizes))
            poles = [Fraction(-1, 3 * size - 4) for size in grouped_sizes if size >= 2]
            grouped_polynomial = secular_coefficients(
                [Fraction(-1, 3 * size - 4) for size in grouped_sizes],
                [Fraction(sizes.count(size) * size, (3 * size - 4) ** 2) for size in grouped_sizes], Fraction(-1))
            assert grouped_polynomial[-1] == 0
            residual = grouped_polynomial[:-1]
            sequence = sturm_sequence(residual)
            endpoints = ['-infinity', *poles]
            intervals = []
            for left, right in zip(endpoints, endpoints[1:]):
                count = sign_variations(sequence, left) - sign_variations(sequence, right)
                assert count == 1
                intervals.append({'left': str(left), 'right': str(right), 'roots': count})
            assert sign_variations(sequence, '-infinity') - sign_variations(sequence, '+infinity') == len(poles)
            negative_spectra.append({'part_sizes': list(sizes),
                                     'residual_polynomial_descending': [str(value) for value in residual],
                                     'root_intervals': intervals})
    assert (distance_cases, inverse_cases, laplacian_polynomials, zero_cofactor_cases) == (2834, 2823, 2818, 16)
    return {'distance_quotient_cases': distance_cases, 'distance_inverse_cases': inverse_cases,
            'laplacian_characteristic_polynomial_cases': laplacian_polynomials,
            'zero_cofactor_distance_cases': zero_cofactor_cases,
            'singular_full_characteristic_polynomials': len(negative_spectra),
            'singular_negative_spectrum_sturm_counts': negative_spectra,
            'checks': ['distance inertia and characteristic polynomial including zero cofactor cases',
                       'distance inverse and Laplacian product identities',
                       'singular distance nullspace and adjugate entry',
                       'parameter bounds and upper equality cases',
                       'Laplacian quotient and full singular characteristic polynomials',
                       'exact Sturm root counts between distinct negative poles']}


def structured_checks():
    examined = defined = singular_inputs = positive_pivots = negative_pivots = 0
    for dimension in range(1, 4):
        for diagonal in product(range(-2, 3), repeat=dimension):
            for vector in product(range(-1, 2), repeat=dimension):
                examined += 1
                pivot = sum(value * coordinate * coordinate for value, coordinate in zip(diagonal, vector))
                if not pivot:
                    continue
                defined += 1
                singular_inputs += 0 in diagonal
                positive_pivots += pivot > 0
                negative_pivots += pivot < 0
                image = [value * coordinate for value, coordinate in zip(diagonal, vector)]
                downdate = [[Fraction(diagonal[row] if row == column else 0) -
                             Fraction(image[row] * image[column], pivot)
                             for column in range(dimension)] for row in range(dimension)]
                expected = [sum(value > 0 for value in diagonal) - (pivot > 0),
                            diagonal.count(0) + 1,
                            sum(value < 0 for value in diagonal) - (pivot < 0)]
                assert inertia(downdate) == expected
                null_vectors = [[int(row == column) for row in range(dimension)]
                                for column, value in enumerate(diagonal) if value == 0] + [list(vector)]
                assert rank(null_vectors) == expected[1]
                assert all(sum(value * coordinate for value, coordinate in zip(values, null_vector)) == 0
                           for values in downdate for null_vector in null_vectors)
    return {'dimension_min': 1, 'dimension_max': 3,
            'diagonal_values': [-2, -1, 0, 1, 2], 'vector_values': [-1, 0, 1],
            'cases_examined': examined, 'nonzero_pivot_cases': defined,
            'singular_input_cases': singular_inputs,
            'positive_pivot_cases': positive_pivots, 'negative_pivot_cases': negative_pivots,
            'checks': ['both structured-downdate inertia triples', 'direct-sum kernel formula']}


def main():
    test_sequence = sturm_sequence(polynomial_from_roots([Fraction(-2), Fraction(-1), Fraction(0), Fraction(1)]))
    test_endpoints = [Fraction(-3), Fraction(-3, 2), Fraction(-1, 2), Fraction(1, 2), Fraction(2)]
    assert all(sign_variations(test_sequence, left) - sign_variations(test_sequence, right) == 1
               for left, right in zip(test_endpoints, test_endpoints[1:]))
    assert sign_variations(test_sequence, '-infinity') - sign_variations(test_sequence, '+infinity') == 4
    assert inertia([[0, 1], [1, 0]]) == [1, 0, 1]
    assert inertia([[0, 0], [0, 0]]) == [0, 2, 0]
    assert inertia([[2, 0], [0, -3]]) == [1, 0, 1]
    assert rank([[0, 1], [0, 0]]) == 1
    examined = defined = excluded = corrected_interval = 0
    singular_examples = []
    for part_count in range(2, 9):
        for sizes in combinations_with_replacement(range(1, 9), part_count):
            if sum(sizes) > 24:
                continue
            examined += 1
            total = sum(sizes)
            nonsingletons = sum(size >= 2 for size in sizes)
            singletons = part_count - nonsingletons
            factors = [3 * size - 4 for size in sizes]
            ratio = sum((Fraction(size, factor) for size, factor in zip(sizes, factors)), Fraction(0))
            if ratio == 0:
                excluded += 1
                continue
            defined += 1
            diagonal, within, between = source_blocks(sizes)
            for row, factor in enumerate(factors):
                assert diagonal[row] == Fraction(1, 4) - Fraction(3, 4 * factor) + 1 / (factor * factor * ratio)
                if sizes[row] >= 2:
                    assert within[row] == -Fraction(3, 4 * factor) + 1 / (factor * factor * ratio)
                for column, other in enumerate(factors):
                    if row != column:
                        assert between[row][column] == 1 / (factor * other * ratio)
            constant_form = [[sizes[row] * diagonal[row] + sizes[row] * (sizes[row] - 1) * within[row]
                              if row == column else sizes[row] * sizes[column] * between[row][column]
                              for column in range(part_count)] for row in range(part_count)]
            constant_inertia = inertia(constant_form)
            expected = ([total - nonsingletons - 1, 1, nonsingletons] if ratio < 0
                        else [total - nonsingletons, 1, nonsingletons - 1])
            assert [constant_inertia[0] + total - part_count,
                    constant_inertia[1], constant_inertia[2]] == expected
            assert rank(shift(constant_form, Fraction(1, 4), sizes)) == part_count
            assert rank(shift(constant_form, Fraction(1), sizes)) == part_count - max(singletons - 1, 0)
            if -1 < ratio < 0:
                corrected_interval += 1
            if ratio != -1:
                continue
            distance, laplacian = full_matrices(sizes)
            assert determinant(distance) == 0
            assert inertia(laplacian) == expected
            assert rank(shift(laplacian, Fraction(1, 4))) == total - (total - part_count)
            assert rank(shift(laplacian, Fraction(1))) == total - (singletons - 1)
            assert all(sum(row) == 0 for row in laplacian)
            cofactor_sum = (-4) ** (total - part_count) * prod(factors) * ratio
            first_minor = [row[1:] for row in laplacian[1:]]
            assert determinant(first_minor) == Fraction((-1) ** (total - 1), cofactor_sum)
            singular_examples.append({'part_sizes': list(sizes), 'inertia': expected,
                                      'quarter_multiplicity': total - part_count,
                                      'one_multiplicity': singletons - 1})
    sizes = (1, 1, 3, 3)
    distance, laplacian = full_matrices(sizes)
    roots = [Fraction(0), Fraction(1)] + [Fraction(1, 4)] * 4 + [Fraction(-1, 5), Fraction(-2)]
    assert charpoly(laplacian) == polynomial_from_roots(roots)
    assert determinant(distance) == 1280
    assert inertia(laplacian) == [5, 1, 2]
    assert rank(shift(laplacian, Fraction(1, 4))) == 4
    print(json.dumps({
        'arithmetic': 'Exact rational arithmetic using Python standard library',
        'range': {'parts_min': 2, 'parts_max': 8, 'part_size_min': 1,
                  'part_size_max': 8, 'vertices_max': 24, 'ordering': 'nondecreasing'},
        'parameter_lists': examined,
        'defined_laplacians_checked': defined,
        'zero_cofactor_lists_excluded': excluded,
        'checks': ['omitted-list entry formulas equal the rank-one representation',
                   'constant-subspace rational congruence inertia',
                   'exact multiplicities of one quarter and one',
                   'singular full-matrix determinant, inertia, eigenspaces, row sums and cofactor'],
        'singular_full_matrix_examples': singular_examples,
        'corrected_interval_examples': corrected_interval,
        'graph_extensions': graph_extension_checks(),
        'structured_downdates': structured_checks(),
        'counterexample': {'part_sizes': list(sizes), 'distance_determinant': 1280,
                          'distance_cofactor_sum': -5120, 'inertia': [5, 1, 2],
                          'characteristic_polynomial': 'x*(x-1)*(x-1/4)^4*(x+1/5)*(x+2)',
                          'descending_coefficients': [str(value) for value in charpoly(laplacian)]},
        'status': 'all exact checks passed'
    }, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
