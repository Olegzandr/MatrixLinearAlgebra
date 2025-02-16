# MatrixLinearAlgebra

Программа на C++ для работы с матрицами и линейными операторами. Реализует вычисление собственных значений/векторов, аффинные преобразования и обработку исключений.

## Особенности
- 🧮 Операции с матрицами: сложение, умножение, вычитание
- 🔍 Вычисление собственных значений и векторов (с использованием библиотеки Eigen)
- 🌀 Геометрические преобразования: поворот, отражение, проекция, масштабирование
- 🛠 Обработка исключений для некорректных операций
- 📁 Сериализация/десериализация матриц через потоки ввода-вывода

## Зависимости
- C++17
- Библиотека [Eigen](https://eigen.tuxfamily.org/) (требуется установка)

## Пример использования
```cpp
#include <iostream>

int main() {
    try {
        BaseMatrix<double> matrix;
        std::cin >> matrix; // Ввод матрицы 2x2
        BaseMatrix<double> result = matrix.eigen_vectors_and_values(matrix);
        std::cout << "Eigen values & vectors:\n" << result;
    } catch (const Exception& e) {
        e.print();
    }
    return 0;
}
```
## Сборка
- Установите Eigen и добавьте путь к заголовочным файлам:
1. git clone https://gitlab.com/libeigen/eigen.git
2. export EIGEN_INCLUDE_PATH=eigen
- Компиляция с поддержкой C++17:
1. g++ -std=c++17 -I$EIGEN_INCLUDE_PATH main.cpp -o matrix_app
- Запустите программу:
1. ./matrix_app

## Обрабатываемые ошибки
- IndexOutOfBoundsException — выход за границы матрицы

- WrongDimensionsException — несовместимые размеры матриц

- NegativeIndexException — отрицательные индексы

- WrongMatrixLinearOperatorException — комплексные решения для вещественных матриц
