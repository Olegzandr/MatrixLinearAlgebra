#include <iostream>
#include <fstream>
#include <Eigen/Eigenvalues>
#include <cstring>
#include <typeinfo>
#include <cmath>


using namespace std;
using namespace Eigen;


class Exception : public std::exception
{
protected:
    //сообщение об ошибке
    char* str;
public:
    Exception(const char* s)
    {
        str = new char[strlen(s) + 1];
        strcpy_s(str, strlen(s) + 1, s);
    }
    Exception(const Exception& e)
    {
        str = new char[strlen(e.str) + 1];
        strcpy_s(str, strlen(e.str) + 1, e.str);
    }
    ~Exception()
    {
        delete[] str;
    }

    virtual void print()
    {
        cout << "\nException: " << str << "\n";
    }
};

class IndexOutOfBoundsException : public Exception
{
protected:
    int index_height;
    int index_width;
public:
    IndexOutOfBoundsException(const char* s, int ind_h, int ind_w) : Exception(s)
    {
        index_height = ind_h;
        index_width = ind_w;
    }

    IndexOutOfBoundsException(const IndexOutOfBoundsException& e) : Exception((Exception&)e)
    {
        index_height = e.index_height;
        index_width = e.index_width;
    }


    void print()
    {
        cout << "\nIndexOutOfBoundsException: " << str << ", non-existent element with indexes: " << index_height << ", " << index_width << "\n";
    }
};

class NegativeIndexException : public IndexOutOfBoundsException
{
public:

    NegativeIndexException(const char* s, int inc_ind) : IndexOutOfBoundsException(s, inc_ind, 0) {};

    NegativeIndexException(const NegativeIndexException& e) : IndexOutOfBoundsException((IndexOutOfBoundsException&)e) {};

    void print()
    {
        cout << "\nNegativeIndexException: " << str << ", negative index: " << index_height << "\n";
    }

};

class TooLargeIndexException : public IndexOutOfBoundsException
{
public:

    TooLargeIndexException(const char* s, int inc_ind) : IndexOutOfBoundsException(s, inc_ind, 0) {};

    TooLargeIndexException(const TooLargeIndexException& e) : IndexOutOfBoundsException((IndexOutOfBoundsException&)e) {};

    void print()
    {
        cout << "\TooLargeIndexException: " << str << ", too large index: " << index_height << "\n";
    }

};

class WrongDimensionsException : public Exception
{
protected:
    int index_height1;
    int index_width1;
    int index_height2;
    int index_width2;
public:
    WrongDimensionsException(const char* s, int index_h1, int index_w1, int index_h2, int index_w2) : Exception(s)
    {
        index_height1 = index_h1;
        index_width1 = index_w1;
        index_height2 = index_h2;
        index_width2 = index_w2;
    }

    WrongDimensionsException(const WrongDimensionsException& e) : Exception((Exception&)e)
    {
        index_height1 = e.index_height1;
        index_width1 = e.index_width1;
        index_height2 = e.index_height2;
        index_width2 = e.index_width2;
    }

    void print()
    {
        cout << "\nWrongDimensionsException: " << str << ", the first matrix has the size: " << index_height1 << "x" << index_width1 << ", the second matrix has the size: " << index_height2 << "x" << index_width2 << "\n";
    }

};

class WrongSizeException : public WrongDimensionsException
{
public:

    WrongSizeException(const char* s, int ind_h, int ind_w) : WrongDimensionsException(s, ind_h, ind_w, 0, 0) {};

    WrongSizeException(const WrongSizeException& e) : WrongDimensionsException((WrongDimensionsException&)e) {};

    void print()
    {
        cout << "\nWrongSizeException: " << str << index_height1 << ", " << index_width1 << "\n";
    }

};

class WrongMatrixLinearOperatorException : public WrongDimensionsException
{
public:

    WrongMatrixLinearOperatorException(const char* s, int ind_h, int ind_w) : WrongDimensionsException(s, ind_h, ind_w, 0, 0) {};

    WrongMatrixLinearOperatorException(const WrongMatrixLinearOperatorException& e) : WrongDimensionsException((WrongDimensionsException&)e) {};

    void print()
    {
        cout << "\WrongMatrixLinearOperatorException: " << str << index_height1 << ", " << index_width1 << "\n";
    }

};

template<class T>
class BaseMatrix;

template<class T>
ostream& operator <<(ostream& ustream, const BaseMatrix<T>& obj)
{
    //запись данных в файл
    if (typeid(ustream).name() == typeid(ofstream).name())
    {
        ustream << obj.height << " " << obj.width << "\n";
        for (int i = 0; i < obj.height; i++)
        {
            for (int j = 0; j < obj.width; j++)
                ustream << obj.ptr[i][j] << " ";
        }
        return ustream;
    }
    //вывод в консоль
    else
    {
        for (int i = 0; i < obj.height; i++)
        {
            for (int j = 0; j < obj.width; j++)
                ustream << obj.ptr[i][j] << " ";
            ustream << "\n";
        }
        return ustream;
    }
}

template<class T>
istream& operator >> (istream& ustream, BaseMatrix<T>& obj)
{
    //считывание данных из файла
    if (typeid(ustream).name() == typeid(ofstream).name())
    {
        for (int i = 0; i < obj.height; i++)
        {
            for (int j = 0; j < obj.width; j++)
                ustream >> obj.ptr[i][j];
        }
    }
    //ввод новых данных о матрице в консоль
    else
    {
        int height_new;
        int width_new;
        ustream >> height_new >> width_new;
        if (height_new == obj.height && width_new == obj.width)
        {
            for (int i = 0; i < obj.height; i++)
            {
                for (int j = 0; j < obj.width; j++)
                    ustream >> obj.ptr[i][j];
            }
            return ustream;
        }
        else
        {
            for (int i = 0; i < obj.height; i++)
                delete[] obj.ptr[i];
            delete[] obj.ptr;
            obj.ptr = NULL;

            if (height_new <= 0 || width_new <= 0)
                throw WrongSizeException("attempt to create a matrix with incorrect sizes: ", height_new, width_new);
            obj.height = height_new;
            obj.width = width_new;
            obj.ptr = new T * [obj.height];
            for (int i = 0; i < height_new; i++)
                obj.ptr[i] = new T[obj.width];

            for (int i = 0; i < obj.height; i++)
            {
                for (int j = 0; j < obj.width; j++)
                    ustream >> obj.ptr[i][j];
            }
            return ustream;
        }
    }
}

template<class T>
class BaseMatrix
{
protected:
    T** ptr;
    int height;
    int width;
public:
    //конструктор
    BaseMatrix(int Height = 2, int Width = 2)
    {
        if (Height <= 0 || Width <= 0)
            throw WrongSizeException("attempt to create a matrix with incorrect sizes: ", Height, Width);
        height = Height;
        width = Width;
        ptr = new T * [height];
        for (int i = 0; i < height; i++)
            ptr[i] = new T[width];
    }

    //конструктор копий
    BaseMatrix(const BaseMatrix& M)
    {
        height = M.height;
        width = M.width;

        ptr = new T * [height];
        for (int i = 0; i < height; i++)
            ptr[i] = new T[width];

        for (int i = 0; i < height; i++)
            for (int j = 0; j < width; j++)
                ptr[i][j] = M.ptr[i][j];
    }

    BaseMatrix operator=(BaseMatrix M)
    {
        if ((width == M.width) and (height == M.height))
        {
            width = M.width;
            height = M.height;

            for (int i = 0; i < height; i++)
                for (int j = 0; j < width; j++)
                    ptr[i][j] = M.ptr[i][j];
        }
        else
        {
            for (int i = 0; i < height; i++)
                delete[] ptr[i];
            delete[] ptr;

            width = M.width;
            height = M.height;

            ptr = new T * [height];
            for (int i = 0; i < height; i++)
            {
                ptr[i] = new T[width];
                for (int j = 0; j < width; j++)
                    ptr[i][j] = M.ptr[i][j];
            }

        }
        return *this;
    }

    //деструктор
    ~BaseMatrix()
    {
        if (ptr != NULL)
        {
            for (int i = 0; i < height; i++)
                delete[] ptr[i];
            delete[] ptr;
            ptr = NULL;
        }
    }

    //вывод
    void print()
    {
        for (int i = 0; i < height; i++)
        {
            for (int j = 0; j < width; j++)
                cout << ptr[i][j] << " ";
            cout << "\n";
        }
    }


    T& operator()(int index1, int index2)
    {
        if (index1 >= height)
            throw TooLargeIndexException("attempt to access an array index that is too large in operator()", index1);
        if (index2 >= width)
            throw TooLargeIndexException("attempt to access an array index that is too large in operator()", index2);
        if (index1 < 0)
            throw NegativeIndexException("attempt to access the negative index of the array in operator()", index1);
        if (index2 < 0)
            throw NegativeIndexException("attempt to access the negative index of the array in operator()", index2);
        return ptr[index1][index2];
    }

    BaseMatrix operator+(BaseMatrix<T> M)
    {
        if (height != M.height || width != M.width)
            throw WrongDimensionsException("attempt to perform an operation between matrices of different sizes in operator+", height, width, M.height, M.width);
        BaseMatrix<T> Res(height, width);
        for (int i = 0; i < height; i++)
            for (int j = 0; j < width; j++)
                Res.ptr[i][j] = ptr[i][j] + M.ptr[i][j];
        return Res;
    }

    BaseMatrix operator-(BaseMatrix<T> M)
    {
        if (height != M.height || width != M.width)
            throw WrongDimensionsException("attempt to perform an operation between matrices of different sizes in operator-", height, width, M.height, M.width);
        BaseMatrix<T> Res(height, width);
        for (int i = 0; i < height; i++)
            for (int j = 0; j < width; j++)
                Res.ptr[i][j] = ptr[i][j] - M.ptr[i][j];
        return Res;
    }

    BaseMatrix operator*(BaseMatrix<T> M)
    {
        if (width != M.height)
            throw WrongDimensionsException("attempt to perform an operation between matrices of incorrect sizes in operator*, in this operation the number of columns in the thirst multiplier is equal to the number of rows in the second", height, width, M.height, M.width);
        BaseMatrix<T> Res(height, M.width);
        int q = width;
        for (int i = 0; i < height; i++)
            for (int j = 0; j < M.width; j++)
            {
                Res.ptr[i][j] = 0;
                for (int k = 0; k < q; k++)
                    Res.ptr[i][j] += (ptr[i][k] * M.ptr[k][j]);
            }

        return Res;
    }

    bool operator==(BaseMatrix<T> M)
    {
        if (height == M.height && width == M.width)
        {
            for (int i = 0; i < height; i++)
                for (int j = 0; j < width; j++)
                {
                    if (ptr[i][j] != M.ptr[i][j])
                        return 0;
                }
        }
        else
            return 0;
        return 1;
    }

    bool operator!=(BaseMatrix<T> M)
    {
        if (height == M.height && width == M.width)
        {
            for (int i = 0; i < height; i++)
                for (int j = 0; j < width; j++)
                {
                    if (ptr[i][j] != M.ptr[i][j])
                        return 1;
                }
        }
        else
            return 1;
        return 0;
    }

    void FastSort(BaseMatrix<T>& obj1, int obj1_first, int obj1_last, BaseMatrix<T>& obj2)
    {
        T temp;
        int f = obj1_first, l = obj1_last;
        T mid = obj1(0, (f + l) / 2);

        do
        {
            while (obj1(0, f) < mid)
                f++;
            while (obj1(0, l) > mid)
                l--;
            if (f <= l)
            {
                //сортировка собственных чисел
                temp = obj1(0, f);
                obj1(0, f) = obj1(0, l);
                obj1(0, l) = temp;

                //сортировка собственных векторов
                for (int i = 0; i <= obj1_last; i++)
                {
                    temp = obj2(i, f);
                    obj2(i, f) = obj2(i, l);
                    obj2(i, l) = temp;
                }

                f++;
                l--;
            }
        } while (f < l);
        {
            //вызов функции рекурсивно,если остались неотсортированные элементы
            if (obj1_first < l)
            {
                FastSort(obj1, obj1_first, l, obj2);
            }
            if (f < obj1_last)
            {
                FastSort(obj1, f, obj1_last, obj2);
            }
        }
    }

    BaseMatrix eigen_vectors_and_values(BaseMatrix<T>& M)
    {
        if (M.height != M.width)
            throw WrongSizeException("attempt to find eigenvalues and eigenvectors from a non-square matrix of a linear operator: ", M.height, M.width);
        //если матрица минимум 3х3
        if (M.width > 2)
        {
            if ((strcmp(typeid(T).name(), "int") == 0) || (strcmp(typeid(T).name(), "double") == 0))
            {
                int row = M.height, col = M.width;
                MatrixXd A(row, col);
                for (int i = 0; i < row; i++)
                {
                    for (int j = 0; j < col; j++)
                    {
                        double d = M(i, j);
                        A(i, j) = d;
                    }
                }
                EigenSolver<MatrixXd> eigensolver(A);
                eigensolver.compute(A, true);
                MatrixXd eigen_values = eigensolver.eigenvalues().real().transpose();
                MatrixXd eigen_vectors = eigensolver.eigenvectors().real();

                BaseMatrix<T> e_values(1, M.width);
                for (int i = 0; i < M.width; i++)
                {
                    double a = eigen_values(i);
                    e_values(0, i) = a;
                }

                BaseMatrix<T> e_vectors(M.width, M.width);
                for (int i = 0; i < M.width; i++)
                {
                    for (int j = 0; j < M.width; j++)
                    {
                        double a = eigen_vectors(i, j);
                        e_vectors(i, j) = a;
                    }
                }

                FastSort(e_values, 0, M.width - 1, e_vectors);

                //матрица с сосбственными числами и векторами
                BaseMatrix<T> e_vectors_and_values(M.width + 1, M.width);
                //заполнение 1 строки(собственными числами)
                for (int i = M.width - 1; i >= 0; i--)
                    e_vectors_and_values(0, M.width - 1 - i) = e_values(0, i);
                //заполнение собственными векторами
                for (int i = 1; i < M.width + 1; i++)
                    for (int j = M.width - 1; j >= 0; j--)
                        e_vectors_and_values(i, M.width - 1 - j) = e_vectors(i - 1, j);

                return e_vectors_and_values;
            }
            else if (strcmp(typeid(T).name(), "float") == 0)
            {
                int row = M.height, col = M.width;
                MatrixXf A(row, col);
                for (int i = 0; i < row; i++)
                {
                    for (int j = 0; j < col; j++)
                    {
                        T d = M(i, j);
                        A(i, j) = d;
                    }
                }
                EigenSolver<MatrixXf> eigensolver(A);
                eigensolver.compute(A, true);
                MatrixXf eigen_values = eigensolver.eigenvalues().real().transpose();
                MatrixXf eigen_vectors = eigensolver.eigenvectors().real();

                BaseMatrix<T> e_values(1, M.width);
                for (int i = 0; i < M.width; i++)
                {
                    T a = eigen_values(i);
                    e_values(0, i) = a;
                }

                BaseMatrix<T> e_vectors(M.width, M.width);
                for (int i = 0; i < M.width; i++)
                {
                    for (int j = 0; j < M.width; j++)
                    {
                        T a = eigen_vectors(i, j);
                        e_vectors(i, j) = a;
                    }
                }

                if (M.width > 1)
                    FastSort(e_values, 0, M.width - 1, e_vectors);

                //матрица с сосбственными числами и векторами
                BaseMatrix<T> e_vectors_and_values(M.width + 1, M.width);
                //заполнение 1 строки(собственными числами)
                for (int i = M.width - 1; i >= 0; i--)
                    e_vectors_and_values(0, M.width - 1 - i) = e_values(0, i);
                //заполнение собственными векторами
                for (int i = 1; i < M.width + 1; i++)
                    for (int j = M.width - 1; j >= 0; j--)
                        e_vectors_and_values(i, M.width - 1 - j) = e_vectors(i - 1, j);

                return e_vectors_and_values;
            }
        }
        if (M.width == 2)
        {
            double det, trace, e_val1, e_val2;
            det = M(0, 0) * M(1, 1) - M(0, 1) * M(1, 0);
            trace = M(0, 0) + M(1, 1);
            e_val1 = trace / 2 + (sqrt(trace * trace / 4 - det));
            e_val2 = trace / 2 - (sqrt(trace * trace / 4 - det));

            double x1, x2, y1, y2;

            //собственные векторы
            if ((strcmp(typeid(T).name(), "float") == 0) || (strcmp(typeid(T).name(), "double") == 0))
            {
                T eps = 0.00000000001;//для сравнения чисел с плавающей точкой

                //проверка на комплексные числа в решении
                if ((fabs(M(0, 0) - M(1, 1)) < eps) && (fabs(M(0, 1) + M(1, 0)) < eps) && ((fabs(M(0, 1)) > eps)) && (fabs(M(1, 0)) > eps))
                    throw WrongMatrixLinearOperatorException("the matrix of a linear operator has complex numbers in the solution: ", M.height, M.width);

                //если матрица диагональная
                if ((fabs(M(0, 0)) > eps) && (M(0, 1) < eps) && (M(1, 0) < eps) && (fabs(M(1, 1)) > eps) && (fabs(M(0, 0) - M(1, 1)) < eps))
                {
                    x1 = 1;
                    y1 = 0;
                    x2 = 0;
                    y2 = 1;
                }

                else
                {
                    if ((fabs(M(0, 1) - 0) < eps) && (fabs(M(0, 0) - e_val1 - 0) < eps))
                    {
                        x1 = M(1, 1) - e_val1;
                        y1 = -M(1, 0);
                    }
                    else
                    {
                        x1 = -M(0, 1);
                        y1 = M(0, 0) - e_val1;
                    }
                    if ((fabs(M(0, 1) - 0) < eps) && (fabs(M(0, 0) - e_val2 - 0) < eps))
                    {
                        x2 = M(1, 1) - e_val2;
                        y2 = -M(1, 0);
                    }
                    else
                    {
                        x2 = -M(0, 1);
                        y2 = M(0, 0) - e_val2;
                    }
                }
            }
            else
            {
                //проверка на комплексные числа в решении
                if ((M(0, 0) - M(1, 1) == 0) && (M(0, 1) + M(1, 0) < 0) && (M(0, 1) > 0) && (M(1, 0) > 0))
                    throw WrongMatrixLinearOperatorException("the matrix of a linear operator has complex numbers in the solution: ", M.height, M.width);
                //если матрица диагональная
                if ((M(0, 0) > 0) && (M(0, 1) < 0) && (M(1, 0) < 0) && (M(1, 1) > 0) && (M(0, 0) - M(1, 1) < 0))
                {
                    x1 = 1;
                    y1 = 0;
                    x2 = 0;
                    y2 = 1;
                } 

                else
                {
                    if ((M(0, 1) == 0) && (M(0, 0) - e_val1 == 0))
                    {
                        x1 = M(1, 1) - e_val1;
                        y1 = -M(1, 0);
                    }
                    else
                    {
                        x1 = -M(0, 1);
                        y1 = M(0, 0) - e_val1;
                    }
                    if ((M(0, 1) == 0) && (M(0, 0) - e_val2 == 0))
                    {
                        x2 = M(1, 1) - e_val2;
                        y2 = -M(1, 0);
                    }
                    else
                    {
                        x2 = -M(0, 1);
                        y2 = M(0, 0) - e_val2;
                    }
                }
            }

            //матрица с сосбственными числами и векторами
            BaseMatrix<T> e_vectors_and_values(3, 2);
            if (e_val1 >= e_val2)
            {
                e_vectors_and_values(0, 0) = e_val1;
                e_vectors_and_values(0, 1) = e_val2;


                e_vectors_and_values(1, 0) = x1;
                e_vectors_and_values(2, 0) = y1;

                e_vectors_and_values(1, 1) = x2;
                e_vectors_and_values(2, 1) = y2;


            }
            else
            {
                e_vectors_and_values(0, 0) = e_val2;
                e_vectors_and_values(0, 1) = e_val1;

                e_vectors_and_values(1, 0) = x2;
                e_vectors_and_values(2, 0) = y2;

                e_vectors_and_values(1, 1) = x1;
                e_vectors_and_values(2, 1) = y1;
            }


            return e_vectors_and_values;
        }

    }

    //поворот
    BaseMatrix Turn(double angle, BaseMatrix<T>& M)
    {
        if (M.height != 2 || M.width != 2)
            throw WrongSizeException("attempt to rotate the matrix on a plane with incorrect sizes: ", M.height, M.width);

        double angle_r = angle * 3.14159265 / 180;//перевод угла в радианы
        //угол неотрицательный(вращение против часовой)
        if (angle >= 0)
        {
            BaseMatrix<T> Res(2, 2);
            Res(0, 0) = cos(angle_r);
            Res(0, 1) = -sin(angle_r);
            Res(1, 0) = sin(angle_r);
            Res(1, 1) = cos(angle_r);
            return (Res * M);
        }
        //(вращение по часовой)
        else
        {
            double angle_r = angle * 3.14159265 / 180;//перевод угла в радианы
            BaseMatrix<T> Res(2, 2);
            Res(0, 0) = cos(angle_r);
            Res(0, 1) = sin(angle_r);
            Res(1, 0) = -sin(angle_r);
            Res(1, 1) = cos(angle_r);
            return (Res * M);
        }
    }

    //отражение(если parametr равен 1 - отражение относительно оси x; если parametr равен 2 - отражение относительно y; 3 - относительно вертикали и горизонтали одновременно;4 - относительно введенной оси(1 параметр) )
    //Матрица К состоит из 2 элементов K(0,0) = a и K(1,0) = b , где y = ax + b
    BaseMatrix Reflection(BaseMatrix<T>& K, int parametr, BaseMatrix<T>& M)
    {
        if (M.height != 2 || M.width != 2)
            throw WrongSizeException("attempt to reflect the matrix on a plane with incorrect sizes: ", M.height, M.width);

        if (parametr == 1)
        {
            BaseMatrix<T> Refl(2, 2);
            Refl(0, 0) = 1;
            Refl(0, 1) = 0;
            Refl(1, 0) = 0;
            Refl(1, 1) = -1;
            return (M * Refl);
        }
        else if (parametr == 2)
        {
            BaseMatrix<T> Refl(2, 2);
            Refl(0, 0) = -1;
            Refl(0, 1) = 0;
            Refl(1, 0) = 0;
            Refl(1, 1) = 1;
            return (M * Refl);
        }
        else if (parametr == 3)
        {
            BaseMatrix<T> Refl(2, 2);
            Refl(0, 0) = -1;
            Refl(0, 1) = 0;
            Refl(1, 0) = 0;
            Refl(1, 1) = -1;
            return (M * Refl);
        }
        else if (parametr == 4)
        {
            if (K.height == 2 && K.width == 1)
            {
                //добавляю единицу к кажому столбцу матрицы линейного оператора
                BaseMatrix<T> New_M(3, 2);
                New_M(0, 0) = M(0, 0);
                New_M(0, 1) = M(0, 1);
                New_M(1, 0) = M(1, 0);
                New_M(1, 1) = M(1, 1);
                New_M(2, 0) = 1;
                New_M(2, 1) = 1;


                //Матрица,которая Перемещает ось так, чтобы она проходила через начало координат
                BaseMatrix<T> Translation_M(3, 3);
                Translation_M(0, 0) = 1;
                Translation_M(0, 1) = 0;
                Translation_M(0, 2) = 0;
                Translation_M(1, 0) = 0;
                Translation_M(1, 1) = 1;
                Translation_M(1, 2) = -K(1, 0);
                Translation_M(2, 0) = 0;
                Translation_M(2, 1) = 0;
                Translation_M(2, 2) = 1;


                //Матрица поворота на угол -tg(a) ( чтобы прямая совпадала с осью х)
                double angle = atan(K(0, 0));//угол в радианах на который нужно повернуть к оси х
                BaseMatrix<T> Turn(3, 3);
                //угол отрицательный(вращение против часовой)
                if (angle <= 0)
                {
                    double angle_pr = -angle;
                    Turn(0, 0) = cos(angle_pr);
                    Turn(0, 1) = -sin(angle_pr);
                    Turn(0, 2) = 0;
                    Turn(1, 0) = sin(angle_pr);
                    Turn(1, 1) = cos(angle_pr);
                    Turn(1, 2) = 0;
                    Turn(2, 0) = 0;
                    Turn(2, 1) = 0;
                    Turn(2, 2) = 1;
                }
                //(вращение по часовой)
                else
                {
                    Turn(0, 0) = cos(angle);
                    Turn(0, 1) = sin(angle);
                    Turn(0, 2) = 0;
                    Turn(1, 0) = -sin(angle);
                    Turn(1, 1) = cos(angle);
                    Turn(1, 2) = 0;
                    Turn(2, 0) = 0;
                    Turn(2, 1) = 0;
                    Turn(2, 2) = 1;
                }


                //Матрица отражения относительно оси х
                BaseMatrix<T> Refl(3, 3);
                Refl(0, 0) = 1;
                Refl(0, 1) = 0;
                Refl(0, 2) = 0;
                Refl(1, 0) = 0;
                Refl(1, 1) = -1;
                Refl(1, 2) = 0;
                Refl(2, 0) = 0;
                Refl(2, 1) = 0;
                Refl(2, 2) = 1;


                //Матрица возвращающая ориентацию прямой на прежнее место(поворачиваем на тот же угол)
                BaseMatrix<T> Turn_Reverse(3, 3);
                //угол отрицательный(вращение по часовой)
                if (angle <= 0)
                {
                    double angle_r = -angle;
                    Turn_Reverse(0, 0) = cos(angle_r);
                    Turn_Reverse(0, 1) = sin(angle_r);
                    Turn_Reverse(0, 2) = 0;
                    Turn_Reverse(1, 0) = -sin(angle_r);
                    Turn_Reverse(1, 1) = cos(angle_r);
                    Turn_Reverse(1, 2) = 0;
                    Turn_Reverse(2, 0) = 0;
                    Turn_Reverse(2, 1) = 0;
                    Turn_Reverse(2, 2) = 1;
                }
                //(вращение против часовой)
                else
                {
                    Turn_Reverse(0, 0) = cos(angle);
                    Turn_Reverse(0, 1) = -sin(angle);
                    Turn_Reverse(0, 2) = 0;
                    Turn_Reverse(1, 0) = sin(angle);
                    Turn_Reverse(1, 1) = cos(angle);
                    Turn_Reverse(1, 2) = 0;
                    Turn_Reverse(2, 0) = 0;
                    Turn_Reverse(2, 1) = 0;
                    Turn_Reverse(2, 2) = 1;
                }


                //Матрица поднимающая прямую обратно на b
                BaseMatrix<T> Translation_M_Reverse(3, 3);
                Translation_M_Reverse(0, 0) = 1;
                Translation_M_Reverse(0, 1) = 0;
                Translation_M_Reverse(0, 2) = 0;
                Translation_M_Reverse(1, 0) = 0;
                Translation_M_Reverse(1, 1) = 1;
                Translation_M_Reverse(1, 2) = K(1, 0);
                Translation_M_Reverse(2, 0) = 0;
                Translation_M_Reverse(2, 1) = 0;
                Translation_M_Reverse(2, 2) = 1;


                BaseMatrix<T> Temp = (Translation_M * Turn * Refl * Turn_Reverse * Translation_M_Reverse * New_M);
                BaseMatrix<T> Res(2, 2);
                Res(0, 0) = Temp(0, 0);
                Res(0, 1) = Temp(0, 1);
                Res(1, 0) = Temp(1, 0);
                Res(1, 1) = Temp(1, 1);
                return Res;

            }
        }
    }

    //проекция(если parametr равен 1 - проекция на ось x; если parametr равен 2 - проекция на ось y)
    BaseMatrix Projection(int parametr, BaseMatrix<T>& M)
    {
        if (M.height != 2 || M.width != 2)
            throw WrongSizeException("attempt to project the matrix on a plane with incorrect sizes: ", M.height, M.width);
        
        if (parametr == 1)
        {
            BaseMatrix<T> Refl(2, 2);
            Refl(0, 0) = 1;
            Refl(0, 1) = 0;
            Refl(1, 0) = 0;
            Refl(1, 1) = 0;
            return (M * Refl);
        }
        else if (parametr == 2)
        {
            BaseMatrix<T> Refl(2, 2);
            Refl(0, 0) = 0;
            Refl(0, 1) = 0;
            Refl(1, 0) = 0;
            Refl(1, 1) = 1;
            return (M * Refl);
        }
    }

    //масштабирование(если parametr = 1 - увеличение в times раз по горизонтали, если parametr = 2 - увеличение в times раз по вертикали,3 - одновременно по вертикали и горизонтали)
    BaseMatrix Scaling(double times, int parametr, BaseMatrix<T>& M)
    {
        if (M.height != 2 || M.width != 2)
            throw WrongSizeException("attempt to scale the matrix on a plane with incorrect sizes: ", M.height, M.width);
        
        if (parametr == 1)
        {
            BaseMatrix<T> Scal(2, 2);
            Scal(0, 0) = times;
            Scal(0, 1) = 0;
            Scal(1, 0) = 0;
            Scal(1, 1) = 1;
            return (M * Scal);
        }
        else if (parametr == 2)
        {
            BaseMatrix<T> Scal(2, 2);
            Scal(0, 0) = 1;
            Scal(0, 1) = 0;
            Scal(1, 0) = 0;
            Scal(1, 1) = times;
            return (M * Scal);
        }
        else if (parametr == 3)
        {
            BaseMatrix<T> Scal(2, 2);
            Scal(0, 0) = times;
            Scal(0, 1) = 0;
            Scal(1, 0) = 0;
            Scal(1, 1) = times;
            return (M * Scal);
        }
    }

    friend ostream& operator<<<T>(ostream&, const BaseMatrix&);
    friend istream& operator>><T>(istream&, BaseMatrix&);
};





int main()
{
    try
    {
        BaseMatrix<double> Q;
        cin >> Q;
        BaseMatrix<double> R = Q.eigen_vectors_and_values(Q);
        cout << "\nResult:\n" << R;
    }

    catch (NegativeIndexException e)
    {
        e.print();
    }
    catch (TooLargeIndexException e)
    {
        e.print();
    }
    catch (IndexOutOfBoundsException e)
    {
        e.print();
    }
    catch (WrongSizeException e)
    {
        e.print();
    }
    catch (WrongMatrixLinearOperatorException e)
    {
        e.print();
    }
    catch (WrongDimensionsException e)
    {
        e.print();
    }
    catch (Exception e)
    {
        e.print();
    }
    catch (exception e)
    {
        cout << "\nAn unknown error has occurred\n";
    }
    

    char c; cin >> c;
    return 0;
}
