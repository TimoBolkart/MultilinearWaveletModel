// diffuse.fs
//
// per-pixel diffuse lighting

varying vec3 N, L0, L1;

void main(void)
{
    const float rhoD = 0.7;

    // output the diffuse color
    float intensity0 = max(0.0, dot(normalize(N), normalize(L0)));
    float intensity1 = max(0.0, dot(normalize(N), normalize(L1)));

    gl_FragColor = gl_Color;
    gl_FragColor.rgb *= rhoD * (intensity0 + intensity1);
}
