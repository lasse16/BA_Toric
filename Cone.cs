using UnityEngine;

internal class Cone
{
    private Vector3 _apex;
    private Vector3 _prefferedVantageAngle;
    //allowed deviation around the preffered vantage angle (degrees) 
    private float _deviationAngle;

    public Cone(Vector3 prefferedVantageAngle, float deviationAngle, Vector3 apex)
    {
        _prefferedVantageAngle = prefferedVantageAngle;
        _deviationAngle = deviationAngle;
        _apex = apex;
    }

    public Cone(Vector3 prefferedVantageAngle, float deviationAngle): this(prefferedVantageAngle,deviationAngle,Vector3.zero)
    {    }

    public Vector3 getMidLineOfCone()
    {
        return _prefferedVantageAngle;
    }

    public float getDiviationAngle()
    {
        return _deviationAngle;
    }

    public bool IsInside(Vector3 testVector)
    {
        return Vector3.Angle(testVector, _prefferedVantageAngle) <= _deviationAngle;
    }

    //TODO check if correct ||check if correct rotation 
    /**
     *Returns the most upward rotated vector around the x axis, 
     * e.g. the *directix* vector rotated by the max deviation angle around the horizontal (x) axis
     * 
     * 
     */
    public Vector3 getBoundaryVectorUp()
    {
        Vector3 res = Quaternion.AngleAxis(- _deviationAngle, Vector3.right) * _prefferedVantageAngle;
        return res;
    }

    public Vector3 getBoundaryVectorDown()
    {
        Vector3 res = Quaternion.AngleAxis(_deviationAngle, Vector3.right) * _prefferedVantageAngle;
        return res;
    }

    public Vector3 getBoundaryVectorLeft()
    {
        Vector3 res = Quaternion.AngleAxis(-_deviationAngle, Vector3.up) * _prefferedVantageAngle;
        return res;
    }

    public Vector3 getBoundaryVectorRight()
    {
        Vector3 res = Quaternion.AngleAxis(_deviationAngle, Vector3.up) * _prefferedVantageAngle;
        return res;
    }




}