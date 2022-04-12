## Test Webapge for Toolbox

This is  matlab toolbox for bla bla

You can download th software package or see examples below.

### Download

You can download the software package [here](https://meyuboglu.github.io/test_webpage/MPC_Poject_2021_Eyuboglu).


For more details see [Basic writing and formatting syntax](https://docs.github.com/en/github/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax).

### Examples

Osicllator. we can put links to code here or direct to another page and explain the example there.
```markdown
Ts = 1;
A = \[1 Ts 0 0; 0 1 0 0; 0 0 1 Ts; 0 0 0 1\];
B = \[-Ts^2/2 0; -Ts 0; Ts^2/2 -Ts^2/2; Ts -Ts\];
n = size(A,1);
m = size(B,2);
C = eye(n);

N = 2;
rho =10

[Link](url) and ![Image](src)
```
Pendulum

Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/meyuboglu/test_webpage/settings/pages). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

### Support or Contact

Having trouble with Pages? Check out our [documentation](https://docs.github.com/categories/github-pages-basics/) or [contact support](https://support.github.com/contact) and we’ll help you sort it out.
