/**
 *  12BCE0478 - ANINDIT KARMAKAR
 *  VIT UNIVERSITY
 *
 *  This program takes as input an Attribute Affinity Matrix, and gives the 
 *  Clustered affinity matrix as its output.
 *  It is upto the user to use the access frequencies from different sites to generate the 
 *  attribute affinity matrix first and then supply that to this program.
 *
 *  The input is taken from the bea.txt file. The first line contains a number N to indicate the
 *  size of the matrix (N x N).
 *  The next N lines contain N space-separated integers. Cumulatively, these integers make up the 
 *  content of the matrix.
 */


#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace std;

int MATRIX_SIZE;

/**
 * The base Matrix class for the required basic matrix operations.
 */
class Matrix {
protected:
    // pointer to hold a 2D array.
    int **_ptr;

    // size of the array
    int _length;

    // temporary 1D array to hold the elements of a column.
    int* _tempColumn;

public:
    /**
     * Constructor to initialize a nxn matrix
     *
     * The content of the matrix starts with index 1.
     * The 0th row contains the column numbers. For example,
     *
     *      A1 A2 A3
     *      10  5  5
     *       6  9  3
     *       2  3  5
     *
     * would be
     *
     *       1  2  3
     *      10  5  5
     *       6  9  3
     *       2  3  5
     */
    Matrix(int n) {
        _length = n;
        _tempColumn = NULL;

        // Initialize the 2D matrix with zeros
        _ptr = (int**) calloc(n+1, sizeof(int*));
        for(int i = 0; i < n+1; i++) {
            _ptr[i] = (int*) calloc(n+1, sizeof(int));
        }

        // Number the colums in the top row. 0th index
        for(int i=1; i<= n; i++) {
            _ptr[0][i] = i;
        }
    }

    /**
     * Inner Proxy class to allow [][] overloading
     *
     * The Matrix class contains a 2D array pointer. Its overloaded [] operator returns
     * an object of Proxy class with contains the 1D array deferenced by the first [] operator.
     *
     * The Proxy class overloads the [] operator to allow the second level of dereferencing
     * with the second [] operator.
     */
    class Proxy {
    private:
        int* _array;

        // Index of the first [] level
        int idx;
    public:
        Proxy(int* arr) : _array(arr) {}

        /**
         * This is the operator overloading function for the second level of [].
         * It returns the final integer value when used like: matrix[5][6]
         */
        int& operator[] (int index) {
            idx = index;
            return _array[index];
        }

        /**
         * This overloads the assignment operator. After the second level [] is dereferenced,
         * it returns a reference to the array element. This assigns the rvalue to the array
         * element.
         */
        void operator= (int rhs) {
            _array[idx] = rhs;
        }
    };

    /**
     * This is the overloading function for the first level of [].
     * It returns a Proxy object containing the 1D array at index;
     */
    Proxy operator[] (int index) {
        return Proxy(_ptr[index]);
    }

    /**
     * This function returns the columns elements of the 2D array at index colIndex.
     * It returns a 1D array containing the column elements for use when copying a column
     * from the attribute affinity matrix to the clustered affinity matrix.
     */
    int* getColumn(int colIndex) {
        free(_tempColumn);
        _tempColumn = (int*) calloc(_length+1, sizeof(int));
        for(int i=0; i<=_length; i++) {
            _tempColumn[i] = _ptr[i][colIndex];
        }

        return _tempColumn;
    }

    /**
     * The name says it all. Prints the contents of the matrix along with the header row.
     */
    void printMatrix() {
        for(int i=0; i<=_length; i++) {
            for(int j=1; j<=_length; j++) {
                cout<<_ptr[i][j]<<"\t";
            }
            cout<<endl;
            if(i==0) {
                for(int k=1;k<=_length; k++)
                    cout<<"--------";
                cout<<endl;
            }
        }
    }
};


/**
 * Derived class for the Attribute Affinity Matrix with some additional functions than the basic
 * Matrix class.
 */
class AttributeMatrix : public Matrix {
public:
    /**
     * Constructor to initilise the matrix.
     */
    AttributeMatrix(int n) : Matrix(n) {}

    /**
     * Function to calculate the bond energy between two columns.
     * Example: bond(A3, A5) would call the function with left=3 and right=5.
     */
    int calculateBond(int left, int right) {
        int sum = 0;
        for(int i = 1; i<=_length; i++) {
            sum = sum + (_ptr[i][left] * _ptr[i][right]);
        }

        return sum;
    }
};


/**
 * Derived class for the Clustered Affinity Matrix with many additional functions than the basic
 * Matrix class.
 */
class ClusteredMatrix : public Matrix {
private:
    // The index of the rightmost column filled in the clustered matrix.
    int _rightmostIndex;

    // When calculating contribution of a placement, we denote is as cont(A1, A2, A5) for example.
    // This means A2 is placed between A1 and A5. Suppose cont(A1, A2, A5) is maximum.

    // contains the column number of the column to the left of the column being placed.
    // (A1 in example)
    int _maxContribLeftIndex;

    // contains the column number of the column that is being placed. (A2 in example)
    int _maxContribMidIndex;

    // contains the column number of the column to the right of the column being placed.
    // (A5 in example)
    int _maxContribRightIndex;

public:

    /**
     * The constructor
     */
    ClusteredMatrix(int n) : Matrix(n) {
        _rightmostIndex = 0;
        _maxContribRightIndex = 0;
        _maxContribMidIndex = 0;
        _maxContribLeftIndex = 0;

        for(int i=1; i<= n; i++) {
            _ptr[0][i] = 0;
        }
    }

    /**
     * Copies the elements of the array in the second parameter to the column denoted by the
     * first parameter.
     */
    void copyToColumn(int colIndex, int* array) {
        for(int i=0; i<=_length; i++) {
            _ptr[i][colIndex] = array[i];
        }
        _rightmostIndex = colIndex;
    }

    /**
     * Getting the index of the rightmost filled column in the clustered affinity matrix.
     */
    int getRightmostIndex() {
        return _rightmostIndex;
    }

    /**
     * Record the three column numbers of the configuration with highest contribution.
     * Example: if cont(A1, A2, A5) has the highest contribution, this function is called with
     * left=1, mid=2, and right=5.
     */
    void recordPlacement(int left, int mid, int right) {
        _maxContribLeftIndex = left;
        _maxContribMidIndex = mid;
        _maxContribRightIndex = right;
    }

    /**
     * Function for placing a column in the clustered matrix from the affinity matrix
     */
    void placeColumnFrom(AttributeMatrix& AA) {

        // Handle placement in leftmost case
        if(_maxContribLeftIndex == 0) {
            // Shift columns to the right
            for(int i = _rightmostIndex+1; i>1; i--) {
                copyToColumn(i, getColumn(i-1));
            }

            // Place the column.
            copyToColumn(1, AA.getColumn(_maxContribMidIndex));

            // Increment the index of the rightmost filled column
            _rightmostIndex++;
            return;
        }

        // For other cases, find the index (start) after which to place the column.
        int start;
        for(start=1; start<=_length; start++) {
            if(_ptr[0][start] == _maxContribLeftIndex)
                break;
        }

        // If start is the rightmost filled index, this is the case when placing in after the
        // rightmost column.
        if(start == _rightmostIndex) {
            _rightmostIndex++;
            copyToColumn(_rightmostIndex, AA.getColumn(_maxContribMidIndex));
        }

        // Handle placement in between two columns.
        // Do shifting.
        for(int i = _rightmostIndex+1; i>start+1; i--) {
            copyToColumn(i, getColumn(i-1));
        }

        // Place the column from the affinity matrix.
        copyToColumn(start+1, AA.getColumn(_maxContribMidIndex));

        // Increment the rightmost filled column index.
        _rightmostIndex++;
    }

    /**
     * Function to make the clustered matrix symmetrical in the final step.
     */
    void makeSymmetrical() {
        int **temp;

        temp = (int**) calloc(_length+1, sizeof(int*));
        for(int i = 0; i < _length+1; i++) {
            temp[i] = (int*) calloc(_length+1, sizeof(int));
        }


        for(int i=1; i<=_length; i++) {
            int row = _ptr[0][i];
            temp[0][i] = row;

            for(int j=1; j<=_length; j++) {
                temp[i][j] = _ptr[row][j];
            }
        }

        // free the previous pointer
        for(int i=0; i<=_length; i++) {
            free(_ptr[i]);
        }
        free(_ptr);

        // Point to new array
        _ptr = temp;
    }
};


/**
 * Function to calculate the bond between two columns
 */
int bond(int left, int right, AttributeMatrix& AA) {
    return AA.calculateBond(left, right);
}


/**
 * Function to calculate the contribution of a certain configuration of columns.
 * Example: to calculate contribution of placement of (A2, A4, A3), this function is called with
 * left = 2, middle = 4, right = 3, and a reference to the Affinity Matrix object.
 */
int calculateContribution(int left, int middle, int right, AttributeMatrix& AA) {
    // if leftmost case
    if(left == 0) {
        return 2*bond(middle, right, AA);
    }

    // if rightmost case
    if(right == middle + 1) {
        return 2*bond(left, middle, AA);
    }

    // if placed in between two columns
    return 2*bond(left, middle, AA) + 2*bond(middle, right, AA) - 2*bond(left, right, AA);
}


/**
 * Function to perform the BEA algorithm.
 */
void doBea(AttributeMatrix& AA, ClusteredMatrix& CA) {

    // Copy the first and second columns from the
    // Attribute Affinity Matrix
    CA.copyToColumn(1, AA.getColumn(1));
    CA.copyToColumn(2, AA.getColumn(2));

    int index = 3;

    while(index <= MATRIX_SIZE) {
        // Declare and Initialize required variables
        int contrib = 0;
        int maxContribution = 0;

        // Calculate placement that gives highest contribution
        for(int i = 1; i < index; i++) {
            contrib = calculateContribution(CA[0][i-1], index, CA[0][i], AA);

            if(contrib >= maxContribution) {
                maxContribution = contrib;
                CA.recordPlacement(CA[0][i-1], index, CA[0][i]);
            }
        }

        contrib = calculateContribution(CA[0][index-1], index, index+1, AA);

        if(contrib >= maxContribution) {
            maxContribution = contrib;
            CA.recordPlacement(CA[0][index-1], index, index+1);
        }

        // Place the column that has the highest contribution
        CA.placeColumnFrom(AA);

        index++;
    }

    CA.makeSymmetrical();
    CA.printMatrix();
}

int main() {
    // Set input stream to the file "bea.txt"
    freopen("bea.txt", "r", stdin);

    // Get the matrix size from the first line of the input file.
    cin>>MATRIX_SIZE;

    // The Attribute Affinity Matrix
    AttributeMatrix AA(MATRIX_SIZE);

    // The Clustered Affinity Matrix
    ClusteredMatrix CA(MATRIX_SIZE);

    // Populate the Attribute Affinity Matrix from the input file.
    // The matrix starts from line 2 onwards.
    for(int i=1; i<=MATRIX_SIZE; i++) {
        for(int j=1; j<=MATRIX_SIZE; j++) {
            int t;
            cin>>t;
            AA[i][j] = t;
        }
    }

    // Do the BEA. Yo!
    doBea(AA, CA);
}
