# ldpc_codec
Low-Density Parity-Check (LDPC) code is a linear error correcting code utilizing for transmitting a message over a noisy transmission channel.

## **Table of Contents**
- [Polar Encoder](#polar-encoder)
- [Polar Decoder](#polar-decoder)
    - [Belief Propagation Summation Product Algorithm](#belief-propagation-summation-product-algorithm)
    - [Min-Summation Product Algorithm](#minsummation-product-algorithm)
    - [Box-Plus Summation Product Algorithm](#box-plus-summation-product-algorithm)

## **Polar Encoder**
## **Polar Decoder**
The commonly procedure for decoding LDPC is the Message Passing Algorithm (MPA), which is built on the Tanner graph with a set of variable nodes (denoted as VP) representing infomation bits and set of check nodes  (denoted by CP) representing parity check bits. The connection between the set of VP nodes and CP nodes follows the parity check matrix (denoted by H) as shown in the figure belows.

![Tanner.png](/image/Tanner.png?raw=true)

MPA procedure is described according to 2 cycles as belows.  
**Initialization Cycle**  
At the initialization cycle, each VP node is assigned a value that equals the value received at the receiver, then each VP node computes the Logarit Likelyhood Ratio (LLR).   
**Loop Cycle**  
Step 1 is **Check-to-Variable messages Processing**. At this step, the CPs compute a Check-to-Variable Message (denoted by CV) to the VPs to which the CP is connected. This calculation is based on information that CPs receive from VPs that are linked to it according to the parity check matrix H.  
Step 2 is **Variable-to-Check messages Processing**. At this step, the VPs compute a Variable-to-Check Message (denoted by VC) to the CPs that the VP is connected. This calculation is based on information that VPs receive from CPs linked to it according to the parity check matrix H.  
Step 3 is **VP Updating**. At this step, the LLR value at each VP is updated based on the information received by the VP messages at that node and the initial LLR value at the initialization step. The result is the algebraic addition of the above values.  
Step 4 is **Hard decision**. At this step, the result obtained after Step 3 is compared with the value 0 to decide whether the transmitted bit is 0 or 1. If the value obtained in Step 3 is greater than or equal to 0, we can say that, the probability that the received bit on the receiver is 0 is higher than the probability that the received bit on the receiver is 1, and decides that the received bit is 0. Otherwise, if the value obtained in Step 3 is less than 0, we can say that, the probability that the bit received on the receiver is 0 is lower than the probability that the received bit on the receiver is 1 and decides that the obtained bit is 1.  
Step 5 is **Syndrome Checking**. The bit sequence obtained after Step 4 is multiplied in a finite field GF(2) by a parity check matrix H. If the gained result is zeros vector, then the bit string obtained in Step 4 is the last decoded bit string after the loop is ended. If the result is a non zeros vector, go back to Step 1 of this cycle.

### **Belief Propagation Summation Product Algorithm**
### **Min-Summation Product Algorithm**
### **Box-Plus Summation Product Algorithm**

## **Reference**
- [1]
- [2]