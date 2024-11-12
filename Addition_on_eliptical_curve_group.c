#include <stdio.h>
#include <stdint.h>

#define NUM_BLOCKS 9
#define BLOCK_BITS 29
#define BLOCK_MASK ((1U << BLOCK_BITS) - 1)  // Mask for 29 bits (0x1FFFFFFF)


void schoolbook_mult(uint64_t *A_blocks, uint64_t *B_blocks, uint64_t *result, int A_size, int B_size);
void Barrett_Reduction(uint64_t *result, int result_size, uint64_t *prime, uint64_t *mu, uint64_t *reduction, int k);
int compare_blocks(uint64_t *a, int a_size, uint64_t *b, int b_size);
void subtract_blocks(uint64_t *a, int a_size, uint64_t *b, int b_size);
void exponential_1(uint64_t *g, int g_size, uint64_t *power, uint64_t *prime, uint64_t *mu, uint64_t *exp_result, int k);
void convert_to_binary(uint64_t *blocks, char *output);


// Schoolbook multiplication of two 256-bit numbers represented as 29-bit blocks
void schoolbook_mult(uint64_t *A_blocks, uint64_t *B_blocks, uint64_t *result, int A_size, int B_size) {
    // Initialize result array with zeros
    for (int i = 0; i < A_size + B_size; i++) {
        result[i] = 0;
    }

    // Perform schoolbook multiplication
    for (int i = 0; i < A_size; i++) {
        for (int j = 0; j < B_size; j++) {
            uint64_t prod = (uint64_t)A_blocks[i] * B_blocks[j];
            result[i + j] += prod & BLOCK_MASK;  // Add lower 29 bits

            // Handle carry
            result[i + j + 1] += prod >> BLOCK_BITS;
        }
    }

    // Propagate carry across result blocks
    for (int i = 0; i < A_size + B_size; i++) {
        uint64_t carry = result[i] >> BLOCK_BITS;
        result[i] &= BLOCK_MASK;  // Keep only the lower 29 bits
        result[i + 1] += carry;
    }
}


void Barrett_Reduction(uint64_t *result, int result_size, uint64_t *prime, uint64_t *mu, uint64_t *reduction, int k) {
    uint64_t q1[k+1];       // To store the top part of result
    uint64_t q2[2 * k + 2];   // Intermediate product q1 * mu
    uint64_t q3[k+1];       // To store q2 divided by prime
    uint64_t q4[2 * k +1];   // Intermediate product q3 * prime
    if(compare_blocks(result, result_size, prime, k) < 0){
      for (int i = 0; i < k ; i++) {
        reduction[i] = result[i];
    }
      return;
    }
    else
    {
	    // Step 1: q1 = top k+1 bits of result (from position k-1 to 2*k-1)
	    for (int i = 0; i < k+1 ; i++) {
		q1[i] = result[i + k - 1];
	    }

	    // Step 2: q2 = q1 * mu (multiply the top part by precomputed mu)
	    schoolbook_mult(q1, mu, q2, k+1, k + 1);  // q2 has k + 1 + k + 1 = 2*k + 2 blocks

	    // Step 3: q3 = top k + 1 bits of q2 (from position k+1 to 2*k+1)
	    for (int i = 0; i < k + 1; i++) {
		q3[i] = q2[i + k + 1];
	    }
	    
	    // Step 4: q4 = q3 * prime (multiply the truncated result by prime)
	    schoolbook_mult(prime, q3, q4, k , k+1);  // q4 has 2*k + 1 blocks
	    
	    // Step 5: reduction = result - q4
	    for (int i = 0; i < result_size; i++) {
		reduction[i] = result[i];
	    }
	    subtract_blocks(reduction, result_size, q4, 2*k+1);
	    
	    // Step 6: Ensure the result is within the prime's range
	    while (compare_blocks(reduction, k, prime, k) >= 0) {
		subtract_blocks(reduction, k, prime, k);
	    }
     }	    	    
}

// Helper function to compare two block arrays (returns 1 if a > b, -1 if a < b, 0 if equal)
int compare_blocks(uint64_t *a, int a_size, uint64_t *b, int b_size) {
    int size = (a_size > b_size) ? a_size : b_size;
    for (int i = size - 1; i >= 0; i--) {
        uint64_t a_val = (i < a_size) ? a[i] : 0;
        uint64_t b_val = (i < b_size) ? b[i] : 0;
        if (a_val > b_val) return 1;
        if (a_val < b_val) return -1;
    }
    return 0;
}


// Helper function to subtract b from a, assuming a >= b
void subtract_blocks(uint64_t *a, int a_size, uint64_t *b, int b_size) {
    int borrow = 0;
    for (int i = 0; i < a_size; i++) {
        uint64_t bi = (i < b_size) ? b[i] : 0;
        if (a[i] < bi + borrow) {
            a[i] = (a[i] + (BLOCK_MASK + 1)) - bi - borrow;
            borrow = 1;
        } else {
            a[i] -= bi + borrow;
            borrow = 0;
        }
    }
}


void convert_to_binary(uint64_t *blocks, char *output) {
    int pos = 0;
    for (int i = 0; i < NUM_BLOCKS; i++) {
        for (int j = BLOCK_BITS - 1; j >= 0; j--) {
            output[pos++] = ((blocks[i] >> j) & 1) + '0';  // Extract each bit and convert to character
        }
    }
    output[pos] = '\0';  // Null-terminate the string
}


// Function for modular exponentiation using Barrett Reduction (right to left)
void exponential_1(uint64_t *g, int g_size, uint64_t *power, uint64_t *prime, uint64_t *mu, uint64_t *exp_result, int k) {
    uint64_t h[18] = {0};  // Fixed size, assuming worst-case for Barrett reduction (max size: prime's size)
    uint64_t temp[36] = {0};  // Temp for multiplication result (2*k blocks)
    uint64_t temp_red[18] = {0};  // Temp for Barrett reduction result (k+1 blocks)
    char binary[NUM_BLOCKS * BLOCK_BITS + 1];
    
    // Convert Large integer into binary
    convert_to_binary(power, binary);
    
    h[0] = 1;  // Initialize h as 1 in modular arithmetic
    for(int j = NUM_BLOCKS * BLOCK_BITS - 1; j >=0; j--){
        if (binary[j] == '1') {
            schoolbook_mult(g, h, temp, g_size, k);
            Barrett_Reduction(temp, g_size + k, prime, mu, temp_red, k);

            // Copy temp_red to h
           // printf("h :  ");
            for (int i = 0; i < k; i++) {
                h[i] = temp_red[i];
               // printf("%ld, ",h[i]);
                temp_red[i] = 0;  // Reset
            }
            //printf("\n\n");
        }

        schoolbook_mult(g, g, temp, g_size, g_size);
        Barrett_Reduction(temp, 2 * g_size, prime, mu, temp_red, k);
        g_size = k;
        // Copy temp_red to g
        //printf("g :  ");
        for (int i = 0; i < k; i++) {
            g[i] = temp_red[i];
            //printf("%ld, ", g[i]);
            temp_red[i] = 0;  // Reset
        }
        //printf("\n\n");
       
    }

    // Copy result to exp_result
    for (int i = 0; i < k; i++) {
        exp_result[i] = h[i];
    }
}

      
      
      
                                                               //     Addition operation on Eliptic Curve Group
                                                               
                                                               
                                                               
      

// Helper function to add b to a, assuming a has enough space to store the result
void add(uint64_t *a, uint64_t *b, uint64_t *result, int size) {
    uint64_t carry = 0;
    for (int i = 0; i < size; i++) {
        uint64_t sum = a[i] + b[i] + carry;
        result[i] = sum & BLOCK_MASK;  // Mask to fit within BLOCK_BITS bits
        carry = sum >> BLOCK_BITS;     // Carry is the overflow beyond BLOCK_BITS
    }
}


// Helper function to subtract b from a, assuming a >= b
void subtract_block(uint64_t *a, uint64_t *b, uint64_t *result, int size) {
    int borrow = 0;
    for (int i = 0; i < size; i++) {
        uint64_t bi = (i < size) ? b[i] : 0;
        if (a[i] < bi + borrow) {
            result[i] = (a[i] + (BLOCK_MASK + 1)) - bi - borrow;
            borrow = 1;
        } else {
            result[i] = a[i] - bi - borrow;
            borrow = 0;
        }
    }
}

// Defined the inverse of an element of eliptic curve group
void inverse(uint64_t *x, uint64_t *inv_res, uint64_t *prime, uint64_t *mu, int k){
   uint64_t x_2[1] = {2};
   uint64_t power[9];
   for(int i=0; i< 9; i++){
      power[i] = prime[i];
   }
   subtract_block(prime, x_2, power, 9); // power = prime - 2
   exponential_1(x, 9, power, prime, mu, inv_res, k); //x^(p-2)
}

// Checks if all elements in both arrays are zero
int check_all_zero(uint64_t *x1, uint64_t *y1, int k) {
    for (int i = 0; i < k; i++) {
        if (x1[i] != 0 || y1[i] != 0) {
            return 0;
        }
    }
    return 1; // If each block of x1 and x2 are zero
}

// Point Addition on Elliptic Curve
void Addition(uint64_t *x1, uint64_t *y1, uint64_t *x2, uint64_t *y2, uint64_t *x3, uint64_t *y3, uint64_t *prime, uint64_t *mu, int k) {
    uint64_t add_result_y[9] = {0};
    uint64_t add_result_x[9] = {0};
    uint64_t sub_result_x[9] = {0};
    uint64_t sub_result_y[9] = {0};
    uint64_t inv_res[9] = {0};
    uint64_t temp[18] = {0};
    uint64_t m[9] = {0};
    uint64_t m_temp[9] = {0};
    uint64_t a_temp[9] = {0,0,0,0,0,0,0,0,-4};

    add(y1, y2, add_result_y, k);   // y1 + y2
    add(x1, x2, add_result_x, k);   // x1 + x2
    subtract_block(x2, x1, sub_result_x, k);  // x2 - x1
    subtract_block(y2, y1, sub_result_y, k);  // y2 - y1

    // Case: Either point is at infinity
    if (check_all_zero(x1, y1, k)) {
        for (int i = 0; i < k; i++) {
            x3[i] = x2[i];
            y3[i] = y2[i];
        }
        return;
    }
    if (check_all_zero(x2, y2, k)) {
        for (int i = 0; i < k; i++) {
            x3[i] = x1[i];
            y3[i] = y1[i];
        }
        return;
    }

    // Case: (y1 + y2 == 0) and (x1 == x2), results in the point at infinity
    if (check_all_zero(add_result_y, add_result_y, k) && check_all_zero(sub_result_x, sub_result_x, k)) {
        for (int i = 0; i < k; i++) {
            x3[i] = 0;
            y3[i] = 0;
        }
        return;
    }

    // Case: x1 != x2
    if (!check_all_zero(sub_result_x, sub_result_x, k)) {
        //compute m
        inverse(sub_result_x, inv_res, prime, mu, k);
        schoolbook_mult(sub_result_y, inv_res, temp, k, k);
        Barrett_Reduction(temp, 2 * k, prime, mu, m, k);
        
        // Compute m^2 - (x2 + x1)
        schoolbook_mult(m, m, temp, k, k);
        Barrett_Reduction(temp, 2 * k, prime, mu, m_temp, k);
        subtract_block(m_temp, add_result_x, x3, k);
        
        // Compute m * (x1 - x3) - y1
        subtract_block(x1, x3, temp, k);
        schoolbook_mult(m, temp, m_temp, k, k);
        Barrett_Reduction(m_temp, 2 * k, prime, mu, temp, k);
        subtract_block(temp, y1, y3, k);
    }
    // case : x1 == x2 and y1 == y2
    else if (check_all_zero(sub_result_x, sub_result_x, k) && check_all_zero(sub_result_y, sub_result_y, k)) {
        add(y1, y1, temp, k); // 2 * y1
        inverse(temp, inv_res, prime, mu, k); // inverse of (2*y1)

        schoolbook_mult(x1, x1, temp, k, k); // x1^2
        Barrett_Reduction(temp, 2 * k, prime, mu, m_temp, k);
        add(m_temp, m_temp, temp, k); // 3 * x1^2
        add(temp, m_temp, temp, k);
        add(temp, a_temp, temp,k); // 3*(x1)^2 + a
        schoolbook_mult(temp, inv_res, m_temp, k, k); //(3*(x1)^2 + a) * inv(2*y1) ---> m
        Barrett_Reduction(m_temp, 2 * k, prime, mu, m, k);
        
        // Compute m^2 - (x2 + x1)
        schoolbook_mult(m, m, temp, k, k);
        Barrett_Reduction(temp, 2 * k, prime, mu, m_temp, k);
        add(x1, x1, temp, k); // 2 * x1
        subtract_block(m_temp, temp, x3, k);
        
        // Compute m * (x1 - x3) - y1
        subtract_block(x1, x3, temp, k);
        schoolbook_mult(m, temp, m_temp, k, k);
        Barrett_Reduction(m_temp, 2 * k, prime, mu, temp, k);
        subtract_block(temp, y1, y3, k);
    }
}      
      


int main() {
    uint64_t prime[9] = {535425013, 174332635, 444665496, 192778653, 388389189, 518147849, 304619691, 363717891, 15281728};
    uint64_t mu[10] = {450887704, 490307913, 387807083, 403879883, 291135210, 307268612, 110539282, 24605042, 70628772, 35};  
    
    // Example curve point (x1, y1) and (x2, y2)
    uint64_t x1[9] = {1,2,2,0,0,0,0,0,0};
    uint64_t y1[9] = {1,2,4,0,0,0,0,0,0};

    uint64_t x2[9] = {1,2,2,0,0,0,0,0,0};
    uint64_t y2[9] = {1,2,4,0,0,0,0,0,0};

    uint64_t x3[9] = {0};
    uint64_t y3[9] = {0};
    int size = 9;

    Addition(x1, y1, x2, y2, x3, y3, prime, mu, size);

    // Print results
    printf("x3: ");
    for (int i = 0; i < size; i++) {
        printf("%ld ", x3[i]);
    }
    printf("\ny3: ");
    for (int i = 0; i < size; i++) {
        printf("%ld ", y3[i]);
    }
    printf("\n");
    
    return 0;
}

