"""
Provides utils function to deal with specific vector and matrix operations.
"""

import numpy as np


def __matrix_multiply(x, y):
    """
    calculates the matrix product for two 1D arrays for the given order.
    :param x: np.array (1D)
    :param y: np.array (1D)
    :return: np.array (2D)  size:(len(x), len(y))
    """
    a = x.reshape(1, len(x))
    b = y.reshape(1, len(y))

    return np.matmul(a.T, b)


def __row_by_vector(x, y, operand):
    """
    Perform selected operation between matrix m-raw and the m-element of
    the vector. Do it for all raw.

    :param x: np.array (2D).
    :param y: np.array (1D).
    :param operand: operand sign
    :return: np.array (2D)
    """

    if operand == '+':
        return x + y.reshape(1, len(y)).T

    elif operand == '-':
        return x - y.reshape(1, len(y)).T

    elif operand == '*':
        return x * y.reshape(1, len(y)).T

    elif operand == '/':
        return x / y.reshape(1, len(y)).T


def __column_by_vector(x, y, operand):
    """
    Perform selected operation between matrix m-column and the m-element of
    the vector. Do it for all columns.

    :param x: np.array (2D).
    :param y: np.array (1D).
    :param operand: operand sign.
    :return: np.array (2D).
    """

    if operand == '+':
        return x + y.reshape(1, len(y))

    elif operand == '-':
        return x - y.reshape(1, len(y))

    elif operand == '*':
        return x * y.reshape(1, len(y))

    elif operand == '/':
        return x / y.reshape(1, len(y))


def __extend_to_matrix(x, dim):
    """
    Creates a matrix with identical rows.
    :param x: np.array (2D).
    :param dim: dimension of matrix. Float.
    :return: np.array (2D)
    """

    return np.array([x, ] * dim).T
